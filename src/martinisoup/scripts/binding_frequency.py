#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from martinisoup.contact import count_contacts, count_contacts_parallel

def index_protein_segments(proteins):
    """
    Build a segment index dictionary with segid, monomer count,
    and a boolean mask for fast per-type lookup.

    Parameters
    ----------
    proteins : MDAnalysis.AtomGroup

    Returns
    -------
    dict
        {protein_name: {"segid": str, "n_monomers": int, "mask": np.ndarray}}
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
        description="Metabolite–protein binding frequency analysis"
    )
    parser.add_argument("topology", type=str, help="Topology file (e.g., .pdb, .tpr)")
    parser.add_argument("trajectory", type=str, help="Trajectory file (e.g., .xtc)")
    parser.add_argument("--protein_sel", type=str, default=None, help="MDAnalysis selection string for proteins")
    parser.add_argument("--metab_sel", type=str, default=None, help="Selection string for metabolites")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Distance cutoff in Å")
    parser.add_argument("--parallel", action="store_true", help="Run in parallel mode")
    parser.add_argument("--n_workers", type=int, default=4, help="Number of workers for parallel mode")
    parser.add_argument("--chunk_size", type=int, default=100, help="Frames per chunk in parallel mode")
    parser.add_argument("--pickle_out", type=str, default="binding.pkl", help="Path to save pickle output")

    args = parser.parse_args()
    command = ' '.join(sys.argv)

    # Load universe
    u = mda.Universe(args.topology, args.trajectory)

    # sort out the selections for the metabolites and proteins
    if not args.metab_sel:
        metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        metabolite_selection = args.metab_sel
    if not args.protein_sel:
        protein_selection = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        protein_selection = args.protein_sel

    proteins = u.select_atoms(protein_selection)
    metabolites = u.select_atoms(metabolite_selection)

    # Build protein segment index with masks
    protein_types = index_protein_segments(proteins)

    # Run analysis
    if args.parallel:
        counts = count_contacts_parallel(
            args.topology,
            args.trajectory,
            protein_selection,
            metabolite_selection,
            protein_types,
            cutoff_radius=args.cutoff,
            n_workers=args.n_workers,
            chunk_size=args.chunk_size
        )
    else:
        counts = count_contacts(
            u,
            proteins,
            metabolites,
            protein_types,
            cutoff_radius=args.cutoff
        )

    # Calculate metabolite concentrations
    metabolite_counts = {}
    for resname in metabolites.residues.resnames:
        metabolite_counts[resname] = metabolite_counts.get(resname, 0) + 1

    protein_results = {}
    system_protein_counts = {}
    for protein_name, protein_counts in counts.items():
        metabolite_hits = {}

        n_protein_monomers = protein_types[protein_name]['n_monomers']
        system_protein_counts[protein_name] = n_protein_monomers

        for metabolite, c in protein_counts.items():
            # Normalize by number of protein monomers and metabolite count
            # This gives: contacts per protein monomer per metabolite
            normalised_count = c / (n_protein_monomers * metabolite_counts[str(metabolite)])
            metabolite_hits[str(metabolite)] = {'count': c,
                                                'normalised_count': normalised_count}
        # sorting the order of the keys may be useful if someone wants to save results as a json or whatever
        protein_results[protein_name] = dict(sorted(metabolite_hits.items()))

    # results for the system.
    results = {'command': command,
               'protein_results': protein_results,
               'n_metabolites': metabolite_counts,
               'n_proteins': system_protein_counts}

    pickle_path = Path(args.pickle_out)
    with open(pickle_path, "wb") as f:
        pickle.dump(results, f)
    print(f"Saved results to {pickle_path}")

if __name__ == "__main__":
    main()
