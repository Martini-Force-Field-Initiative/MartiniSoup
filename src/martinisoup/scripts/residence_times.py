#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np

from MDAnalysis import Universe
from martinisoup.residence_tracker import BindingEventTracker
from MDAnalysis import transformations

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
    parser.add_argument("--cutoff", type=float, default=5.0,
                        help="Distance cutoff in Angstrom")
    parser.add_argument("--start", type=int, default=0,
                        help="Start frame")
    parser.add_argument("--stop", type=int, default=None,
                        help="Stop frame")
    parser.add_argument("--step", type=int, default=1,
                        help="Frame step size")
    parser.add_argument("--output", default="lifetimes.pkl",
                        help="Output pickle file for molecule-level results")
    parser.add_argument("--summary", action="store_true",
                        help="Store unique residence times and their counts instead of "
                             "the full duration list, producing a smaller output file")

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

    # Load universe
    print("Loading universe...")
    u = Universe(args.topology, args.trajectory)
    # make sure all atoms are in the box
    ag = u.atoms
    transform = transformations.wrap(ag)
    u.trajectory.add_transformations(transform)

    metabolites = u.select_atoms(metabolite_selection)
    proteins = u.select_atoms(protein_selection)

    # Run trajectory analysis
    tracker = BindingEventTracker(u, metabolites, proteins,
                                  cutoff=args.cutoff,
                                  start=args.start,
                                  stop=args.stop,
                                  step=args.step)
    residences = tracker.track()

    # Save full results
    pickle_path = Path(args.output)
    with open(pickle_path, "wb") as f:
        pickle.dump({"command": command, "residences": residences}, f)
    print(f"Saved results to {pickle_path}")

    if args.summary:
        summary = {
            moltype: np.unique(durations, return_counts=True)
            for moltype, durations in residences.items()
        }
        summary_path = pickle_path.with_stem(pickle_path.stem + "_summary")
        with open(summary_path, "wb") as f:
            pickle.dump({"command": command, "residences": summary}, f)
        print(f"Saved summary to {summary_path}")

if __name__ == "__main__":
    main()
