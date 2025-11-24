"""
Command-line interface for metabolite–protein contact analysis
using binding_analysis.contact.
"""

import argparse
import pickle
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from contact import run_contact_analysis, parallel_contact_analysis

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
            chunk_size=args.chunk_size
        )
    else:
        counts = run_contact_analysis(
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
    for protein_name, protein_counts in counts.items():
        metabolite_hits = {}
        n_protein_monomers = protein_types[protein_name]['n_monomers']
        for metabolite, c in protein_counts.items():
            # Normalize by number of protein monomers and metabolite count
            # This gives: contacts per protein monomer per metabolite
            normalised_count = c / (n_protein_monomers * metabolite_counts[str(metabolite)])
            metabolite_hits[str(metabolite)] = {'count': c,
                                                'normalised_count': normalised_count}
        protein_results[protein_name] = {'n_monomers': n_protein_monomers,
                                         'metabolite_counts': metabolite_hits
                                         }

    # results for the system.
    results = {'protein_results': protein_results,
               'metabolite_counts': metabolite_counts}

    # Save pickle if requested
    if args.pickle_out:
        pickle_path = Path(args.pickle_out)
        with open(pickle_path, "wb") as f:
            pickle.dump(results, f)
        print(f"Saved results to {pickle_path}")

if __name__ == "__main__":
    main()
