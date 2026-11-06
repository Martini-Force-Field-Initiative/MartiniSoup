#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

import MDAnalysis as mda

from martinisoup.protein_clustering import (
    analyse_trajectory,
    analyse_trajectory_parallel,
    build_nodes,
)
from martinisoup.units import Q_, UnitSystem, quantity_arg, unit_string


def main():
    parser = argparse.ArgumentParser(
        description="Build per-frame protein contact networks from an MD trajectory."
    )
    parser.add_argument("topology", type=str, help="Topology file (e.g., .tpr, .gro)")
    parser.add_argument("trajectory", type=str, help="Trajectory file (e.g., .xtc, .trr)")
    parser.add_argument("--protein-selection", default=None,
                        help="MDAnalysis selection string for proteins")
    parser.add_argument("--r-max", type=quantity_arg('angstrom'), default=Q_(6.0, 'angstrom'),
                        help="Contact distance cutoff (default: 6 Å). A bare number is in Å; "
                             "a value with units, e.g. '0.6 nm', is converted.")
    parser.add_argument("--start", type=int, default=0,
                        help="Start frame (default: 0)")
    parser.add_argument("--stop", type=int, default=None,
                        help="Stop frame, exclusive (default: end of trajectory)")
    parser.add_argument("--step", type=int, default=1,
                        help="Frame stride (default: 1)")
    parser.add_argument("--output", type=str, default="protein_contacts.pkl",
                        help="Output pickle file (default: protein_contacts.pkl)")
    parser.add_argument("--parallel", action="store_true",
                        help="Run analysis in parallel")
    parser.add_argument("--n_workers", type=int, default=4,
                        help="Number of worker processes for parallel mode (default: 4)")
    parser.add_argument("--chunk_size", type=int, default=100,
                        help="Approximate frames per part in parallel mode (default: 100)")

    args = parser.parse_args()
    command = ' '.join(sys.argv)

    if not args.protein_selection:
        protein_selection = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        protein_selection = args.protein_selection

    if args.parallel:
        frames = analyse_trajectory_parallel(
            args.topology,
            args.trajectory,
            protein_selection,
            r_max=args.r_max,
            start=args.start,
            stop=args.stop,
            step=args.step,
            n_workers=args.n_workers,
            chunk_size=args.chunk_size,
        )
        # Rebuild nodes from topology for the top-level output metadata
        u = mda.Universe(args.topology, args.trajectory)
        nodes = build_nodes(u.select_atoms(protein_selection))
    else:
        u = mda.Universe(args.topology, args.trajectory)
        atoms = u.select_atoms(protein_selection)
        nodes = build_nodes(atoms)
        frames = analyse_trajectory(
            u,
            atoms,
            nodes,
            r_max=args.r_max,
            start=args.start,
            stop=args.stop,
            step=args.step,
        )

    unit_system = UnitSystem.from_universe(u)
    output_path = Path(args.output)
    with open(output_path, "wb") as f:
        pickle.dump({
            'command': command,
            'nodes': nodes,
            'r_max': unit_system.length_magnitude(args.r_max),
            'frames': frames,
            'units': {
                'r_max': unit_string(unit_system.length_unit),
                'frames.frame': 'dimensionless',
                'frames.time': unit_string(unit_system.time_unit),
            },
        }, f)
    print(f"Saved results to {output_path}")


if __name__ == "__main__":
    main()
