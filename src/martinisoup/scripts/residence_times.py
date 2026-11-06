#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from martinisoup.residence_tracker import track_serial, track_parallel
from martinisoup.units import Q_, UnitSystem, quantity_arg, unit_string

def main():
    parser = argparse.ArgumentParser(
        description="Analyze metabolite-protein binding event residence times from an MD trajectory."
    )
    parser.add_argument("topology", help="Topology file (e.g., PDB, PSF, GRO)")
    parser.add_argument("trajectory", help="Trajectory file (e.g., DCD, XTC, TRR)")
    parser.add_argument("--metabolite-selection", default=None,
                        help="MDAnalysis atom selection string for metabolites")
    parser.add_argument("--protein-selection", default=None,
                        help="MDAnalysis atom selection string for proteins")
    parser.add_argument("--cutoff", type=quantity_arg('angstrom'), default=Q_(5.0, 'angstrom'),
                        help="Distance cutoff (default: 5 Å). A bare number is in Å; "
                             "a value with units, e.g. '0.5 nm', is converted.")
    parser.add_argument("--start", type=int, default=0,
                        help="Start frame")
    parser.add_argument("--stop", type=int, default=None,
                        help="Stop frame")
    parser.add_argument("--step", type=int, default=1,
                        help="Frame step size")
    parser.add_argument("--output", default="lifetimes.pkl",
                        help="Output pickle file for molecule-level results")
    parser.add_argument("--frames", action="store_true",
                        help="Measure durations as numbers of trajectory frames instead "
                             "of time (default: time, in ps). The time between frames is "
                             "saved as 'dt', so durations can still be converted to time "
                             "(residence-fitter does this automatically).")
    parser.add_argument("--summary", default=None, metavar="FILE",
                        help="Write a summary pickle (unique residence times and their "
                             "counts) to FILE")
    parser.add_argument("--parallel", action="store_true",
                        help="Run contact detection in parallel (Phase 1), then run "
                             "the state machine sequentially (Phase 2). Produces "
                             "identical results to the serial mode.")
    parser.add_argument("--n_workers", type=int, default=4,
                        help="Number of worker processes for parallel mode")
    parser.add_argument("--chunk_size", type=int, default=100,
                        help="Approximate frames per part in parallel mode")

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

    # Run trajectory analysis
    if args.parallel:
        residences = track_parallel(
            args.topology,
            args.trajectory,
            metabolite_selection,
            protein_selection,
            cutoff=args.cutoff,
            start=args.start,
            stop=args.stop,
            step=args.step,
            n_workers=args.n_workers,
            chunk_size=args.chunk_size,
            use_time=not args.frames,
        )
    else:
        residences = track_serial(
            args.topology,
            args.trajectory,
            metabolite_selection,
            protein_selection,
            cutoff=args.cutoff,
            start=args.start,
            stop=args.stop,
            step=args.step,
            use_time=not args.frames,
        )

    # Units of the data as read (track_* read the files with default settings)
    u = mda.Universe(args.topology, args.trajectory)
    unit_system = UnitSystem.from_universe(u)
    time_unit = unit_string(unit_system.time_unit)
    duration_unit = 'dimensionless' if args.frames else time_unit
    cutoff = unit_system.length_magnitude(args.cutoff)
    length_unit = unit_string(unit_system.length_unit)
    common = {
        "stamp": "frame" if args.frames else "time",
        "dt": u.trajectory.dt,  # time between frames
        "cutoff": cutoff,
    }

    # Save full results
    pickle_path = Path(args.output)
    with open(pickle_path, "wb") as f:
        pickle.dump({
            "command": command,
            "residences": residences,
            **common,
            "units": {"residences": duration_unit, "dt": time_unit, "cutoff": length_unit},
        }, f)
    print(f"Saved results to {pickle_path}")

    if args.summary:
        summary = {
            moltype: np.unique(durations, return_counts=True)
            for moltype, durations in residences.items()
        }
        summary_path = Path(args.summary)
        with open(summary_path, "wb") as f:
            pickle.dump({
                "command": command,
                "residences": summary,  # {moltype: (values, counts)}
                **common,
                "units": {
                    "residences.values": duration_unit,
                    "residences.counts": "dimensionless",
                    "dt": time_unit,
                    "cutoff": length_unit,
                },
            }, f)
        print(f"Saved summary to {summary_path}")

if __name__ == "__main__":
    main()
