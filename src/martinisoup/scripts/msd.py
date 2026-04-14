#!/usr/bin/env python3

import argparse
import pickle
import sys
from pathlib import Path

from martinisoup.msd import compute_msd


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

    metab_sel = args.metab_sel or (
        'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE '
        'GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    )

    print("Loading universe...")
    results = compute_msd(
        args.topology,
        args.trajectory,
        metab_sel,
        start=args.start,
        stop=args.stop,
        fft=args.fft,
    )
    results['command'] = command

    output_path = Path(args.output)
    with open(output_path, 'wb') as f:
        pickle.dump(results, f)
    print(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
