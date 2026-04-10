#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

import MDAnalysis as mda

from martinisoup.clustering import analyse_trajectory, analyse_trajectory_parallel


def main():
    parser = argparse.ArgumentParser(
        description="Classify metabolites as clustered, soluble, or protein-adsorbed across an MD trajectory."
    )
    parser.add_argument("topology", type=str, help="Topology file (e.g., .tpr, .gro)")
    parser.add_argument("trajectory", type=str, help="Trajectory file (e.g., .xtc, .trr)")
    parser.add_argument("--metabolite-selection", default=None,
                        help="MDAnalysis selection string for metabolites")
    parser.add_argument("--protein-selection", default=None,
                        help="MDAnalysis selection string for proteins")
    parser.add_argument("--r-max", type=float, default=5.0,
                        help="Cluster neighbour cutoff in Angstrom (default: 5.0)")
    parser.add_argument("--contact-cutoff", type=float, default=5.0,
                        help="Protein contact cutoff in Angstrom (default: 5.0)")
    parser.add_argument("--start", type=int, default=0,
                        help="Start frame (default: 0)")
    parser.add_argument("--stop", type=int, default=None,
                        help="Stop frame, exclusive (default: end of trajectory)")
    parser.add_argument("--step", type=int, default=1,
                        help="Frame stride (default: 1)")
    parser.add_argument("--output", type=str, default="cluster_states.pkl",
                        help="Output pickle file (default: cluster_states.pkl)")
    parser.add_argument("--parallel", action="store_true",
                        help="Run analysis in parallel")
    parser.add_argument("--n_workers", type=int, default=4,
                        help="Number of worker processes for parallel mode (default: 4)")
    parser.add_argument("--chunk_size", type=int, default=100,
                        help="Frames per chunk in parallel mode (default: 100)")

    args = parser.parse_args()
    command = ' '.join(sys.argv)

    if not args.metabolite_selection:
        metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        metabolite_selection = args.metabolite_selection

    if not args.protein_selection:
        protein_selection = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        protein_selection = args.protein_selection

    if args.parallel:
        results = analyse_trajectory_parallel(
            args.topology,
            args.trajectory,
            metabolite_selection,
            protein_selection,
            r_max=args.r_max,
            contact_cutoff=args.contact_cutoff,
            start=args.start,
            stop=args.stop,
            step=args.step,
            n_workers=args.n_workers,
            chunk_size=args.chunk_size,
        )
    else:
        u = mda.Universe(args.topology, args.trajectory)
        metabolites = u.select_atoms(metabolite_selection)
        proteins = u.select_atoms(protein_selection)

        results = analyse_trajectory(
            u,
            metabolites,
            proteins,
            r_max=args.r_max,
            contact_cutoff=args.contact_cutoff,
            start=args.start,
            stop=args.stop,
            step=args.step,
        )

    output_path = Path(args.output)
    with open(output_path, "wb") as f:
        pickle.dump({"command": command, "results": results}, f)
    print(f"Saved results to {output_path}")


if __name__ == "__main__":
    main()
