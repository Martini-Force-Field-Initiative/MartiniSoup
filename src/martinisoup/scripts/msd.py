#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd_analysis
import numpy as np
from MDAnalysis.transformations import NoJump
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser(
        description="Mean squared displacement analysis per metabolite type."
    )
    parser.add_argument("topology", help="Topology file (e.g., .tpr, .gro)")
    parser.add_argument("trajectory", help="Trajectory file (e.g., .xtc, .trr)")
    parser.add_argument("--metab_sel", type=str, default=None,
                        help="MDAnalysis selection string for metabolites")
    parser.add_argument("--start", type=int, default=None,
                        help="Start frame (default: first frame of second half of trajectory)")
    parser.add_argument("--stop", type=int, default=None,
                        help="Stop frame (default: start + 100)")
    parser.add_argument("--fft", action="store_true", default=True,
                        help="Use FFT algorithm for MSD calculation (default: True)")
    parser.add_argument("--output", type=str, default="msd.pkl",
                        help="Output pickle file (default: msd.pkl)")

    args = parser.parse_args()
    command = ' '.join(sys.argv)

    if not args.metab_sel:
        metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        metabolite_selection = args.metab_sel

    print("Loading universe...")
    u = mda.Universe(args.topology, args.trajectory)
    u.trajectory.add_transformations(NoJump())

    # Default frame range: second half of trajectory, up to 100 frames
    start = args.start if args.start is not None else int(len(u.trajectory) / 2)
    stop = args.stop if args.stop is not None else start + 100

    print(f"Running MSD on frames {start}–{stop}...")
    MSD = msd_analysis.EinsteinMSD(u, select=metabolite_selection, fft=args.fft)
    MSD.run(start=start, stop=stop, verbose=True)

    msd_data = {key: value for key, value in MSD.results.items()}
    msd_array = msd_data['msds_by_particle']

    metabolites = u.select_atoms(metabolite_selection)
    unique_resnames = list(set(metabolites.resnames))

    residue_average = np.zeros((len(unique_resnames), msd_array.shape[0]))
    residue_std = np.zeros((len(unique_resnames), msd_array.shape[0]))

    for i, resname in enumerate(tqdm(unique_resnames, desc="Averaging per molecule type")):
        mask = metabolites.resnames == resname
        residue_average[i] = np.mean(msd_array.T[mask], axis=0)
        residue_std[i] = np.std(msd_array.T[mask], axis=0)

    results = {
        'command': command,
        'resnames': unique_resnames,
        'residue_timeseries': residue_average,
        'residue_std': residue_std,
        'time': MSD.times,
        'lagtimes': msd_data['timeseries'],
        'dimensions': MSD.dim_fac,
        'dt': u.trajectory.dt,
    }

    output_path = Path(args.output)
    with open(output_path, 'wb') as f:
        pickle.dump(results, f)
    print(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
