#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

import MDAnalysis as mda

from martinisoup.protein_rdf import compute_rdf, compute_rdf_parallel


def main():
    parser = argparse.ArgumentParser(
        description="Compute protein centre-of-mass RDF from an MD trajectory."
    )
    parser.add_argument("topology",   type=str, help="Topology file (e.g., .tpr, .gro)")
    parser.add_argument("trajectory", type=str, help="Trajectory file (e.g., .xtc, .trr)")
    parser.add_argument("--protein-selection", default=None,
                        help="MDAnalysis selection string for protein atoms")
    parser.add_argument("--r-max", type=float, default=100.0,
                        help="Upper distance limit in Angstrom (default: 100.0)")
    parser.add_argument("--n-bins", type=int, default=200,
                        help="Number of histogram bins (default: 200)")
    parser.add_argument("--start", type=int, default=0,
                        help="Start frame (default: 0)")
    parser.add_argument("--stop", type=int, default=None,
                        help="Stop frame, exclusive (default: end of trajectory)")
    parser.add_argument("--step", type=int, default=1,
                        help="Frame stride (default: 1)")
    parser.add_argument("--output", type=str, default="protein_rdf.pkl",
                        help="Output pickle file (default: protein_rdf.pkl)")
    parser.add_argument("--parallel", action="store_true",
                        help="Run analysis in parallel")
    parser.add_argument("--n_workers", type=int, default=4,
                        help="Number of worker processes for parallel mode (default: 4)")
    parser.add_argument("--chunk_size", type=int, default=100,
                        help="Frames per chunk in parallel mode (default: 100)")

    args    = parser.parse_args()
    command = ' '.join(sys.argv)

    if args.protein_selection:
        prot_sel = args.protein_selection
    else:
        prot_sel = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'

    if args.parallel:
        result = compute_rdf_parallel(
            args.topology,
            args.trajectory,
            prot_sel,
            r_max=args.r_max,
            n_bins=args.n_bins,
            start=args.start,
            stop=args.stop,
            step=args.step,
            n_workers=args.n_workers,
            chunk_size=args.chunk_size,
        )
    else:
        u     = mda.Universe(args.topology, args.trajectory)
        atoms = u.select_atoms(prot_sel)
        result = compute_rdf(
            u,
            atoms,
            r_max=args.r_max,
            n_bins=args.n_bins,
            start=args.start,
            stop=args.stop,
            step=args.step,
        )

    output_path = Path(args.output)
    with open(output_path, 'wb') as f:
        pickle.dump({'command': command, **result}, f)
    print(f"Saved results to {output_path}")


if __name__ == "__main__":
    main()
