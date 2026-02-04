#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from cytosol_analysis.contact import run_contact_analysis, parallel_contact_analysis

def build_protein_types(proteins):
    """
    Build protein_types dictionary with:
        - segid
        - n_monomers
        - boolean mask for fast lookup
    """
    protein_types = {}
    for seg in proteins.segments:
        # Protein name from segid
        parts = seg.segid.split('_')
        protein_name = '_'.join(parts[2:])

        # Count monomers
        n_monomers = len(np.unique(seg.atoms.molnums))

        # Boolean mask for this protein
        mask = proteins.segids == seg.segid

        protein_types[protein_name] = {
            "segid": seg.segid,
            "n_monomers": n_monomers,
            "mask": mask
        }
    return protein_types


def main():
    parser = argparse.ArgumentParser(
        description="Metabolite–protein contact analysis"
    )
    parser.add_argument("topology", type=str, help="Topology file (e.g., .pdb, .tpr)")
    parser.add_argument("trajectory", type=str, help="Trajectory file (e.g., .xtc)")
    parser.add_argument("--protein_sel", type=str, default=None, help="MDAnalysis selection string for proteins")
    parser.add_argument("--metab_sel", type=str, default=None, help="Selection string for metabolites")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Distance cutoff in Å")
    parser.add_argument("--residues", action="store_true", help="Track metabolites at residue level")
    parser.add_argument("--parallel", action="store_true", help="Run in parallel mode")
    parser.add_argument("--n_workers", type=int, default=4, help="Number of workers for parallel mode")
    parser.add_argument("--chunk_size", type=int, default=100, help="Frames per chunk in parallel mode")
    parser.add_argument("--pickle_out", type=str, default=None, help="Path to save pickle output")

    args = parser.parse_args()

    # Load universe
    u = mda.Universe(args.topology, args.trajectory)

    # sort out the selections for the metabolites and proteins
    if not args.metab_sel:
        metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        metabolite_selection = 'resname ATP'
    if not args.protein_sel:
        protein_selection = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        protein_selection = 'protein'

    proteins = u.select_atoms(protein_selection)
    metabolites = u.select_atoms(metabolite_selection)

    # Build protein_types dict with masks
    protein_types = build_protein_types(proteins)

    # Run analysis
    if args.parallel:
        counts = parallel_contact_analysis(
            args.topology,
            args.trajectory,
            protein_selection,
            metabolite_selection,
            protein_types,
            cutoff_radius=args.cutoff,
            n_workers=args.n_workers,
            chunk_size=args.chunk_size,
            track_residues=args.residues
        )
    else:
        counts = run_contact_analysis(
            u,
            proteins,
            metabolites,
            protein_types,
            cutoff_radius=args.cutoff,
            track_residues=args.residues
        )
    # ------------------------------------------------------------
    # Calculate metabolite concentrations (number of molecules)
    # ------------------------------------------------------------
    metabolite_counts = {}
    for resname in metabolites.residues.resnames:
        metabolite_counts[resname] = metabolite_counts.get(resname, 0) + 1

    # ------------------------------------------------------------
    # Detect result structure
    # ------------------------------------------------------------
    if "atom_level" in counts:
        atom_counts = counts["atom_level"]
        residue_counts = counts.get("residue_level", None)
    else:
        atom_counts = counts
        residue_counts = None

    # ------------------------------------------------------------
    # Protein-level aggregation (unchanged)
    # ------------------------------------------------------------
    protein_results = {}
    system_protein_counts = {}

    for protein_name, protein_counts in atom_counts.items():

        metabolite_hits = {}

        n_protein_monomers = protein_types[protein_name]["n_monomers"]
        system_protein_counts[protein_name] = n_protein_monomers

        for metabolite, c in protein_counts.items():
            n_met = metabolite_counts.get(str(metabolite), 1)
            normalised_count = c / (n_protein_monomers * n_met)

            metabolite_hits[str(metabolite)] = {
                "count": int(c),
                "normalised_count": float(normalised_count),
            }

        protein_results[protein_name] = dict(sorted(metabolite_hits.items()))

    # ------------------------------------------------------------
    # Residue-level aggregation with monomer correction
    # ------------------------------------------------------------
    residue_results = None

    if residue_counts is not None:
        residue_results = {}

        for protein_name, met_dict in residue_counts.items():
            # ------------------------------------------------------------
            # Build mapping: global resid -> (resname, local_resid)
            # ------------------------------------------------------------
            pdata = protein_types[protein_name]
            mask = pdata["mask"]
            prot_atoms = proteins[mask]

            # Work at residue level
            prot_residues = prot_atoms.residues

            # Group residues by molnum
            molnums = prot_residues.molnums
            unique_molnums = np.unique(molnums)

            # Use first monomer as reference
            ref_mol = unique_molnums[0]
            ref_residues = prot_residues[molnums == ref_mol]

            # Assign local residue indices (1-based, in order)
            ref_resid_to_local = {}
            for local_idx, res in enumerate(ref_residues, start=1):
                ref_resid_to_local[res.resid] = (res.resname, local_idx)

            # Build full mapping for all monomers
            resid_to_local = {}

            for mol in unique_molnums:
                mol_residues = prot_residues[molnums == mol]

                if len(mol_residues) != len(ref_residues):
                    raise ValueError(
                        f"Protein {protein_name}: monomer {mol} has "
                        f"{len(mol_residues)} residues, expected {len(ref_residues)}"
                    )

                for ref_res, res in zip(ref_residues, mol_residues):
                    resid_to_local[res.resid] = (ref_res.resname, ref_resid_to_local[ref_res.resid][1])

            protein_residue_data = {}

            for metabolite, res_dict in met_dict.items():
                metabolite_residue_hits = {}

                n_met = metabolite_counts.get(str(metabolite), 1)
                n_protein_monomers = pdata["n_monomers"]

                for (prot_resname, global_resid), c in res_dict.items():

                    # Collapse across monomers
                    if global_resid not in resid_to_local:
                        continue

                    _, local_resid = resid_to_local[global_resid]

                    key = (prot_resname, local_resid)

                    normalised_count = c / (n_protein_monomers * n_met)

                    if key not in metabolite_residue_hits:
                        metabolite_residue_hits[key] = {
                            "count": 0,
                            "normalised_count": 0.0,
                        }

                    metabolite_residue_hits[key]["count"] += int(c)
                    metabolite_residue_hits[key]["normalised_count"] += float(
                        normalised_count
                    )

                protein_residue_data[str(metabolite)] = metabolite_residue_hits

            residue_results[protein_name] = protein_residue_data

    # ------------------------------------------------------------
    # Final results dictionary
    # ------------------------------------------------------------
    results = {
        "protein_results": protein_results,
        "n_metabolites": metabolite_counts,
        "n_proteins": system_protein_counts,
    }

    if residue_results is not None:
        results["residue_results"] = residue_results

    results["command_used"] = str(sys.argv[1:])

    # ------------------------------------------------------------
    # Save pickle if requested
    # ------------------------------------------------------------
    if args.pickle_out:
        pickle_path = Path(args.pickle_out)
        with open(pickle_path, "wb") as f:
            pickle.dump(results, f)
        print(f"Saved results to {pickle_path}")


if __name__ == "__main__":
    main()
