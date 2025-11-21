"""
Command-line interface for binding_analysis.
"""

import argparse
import pickle
from MDAnalysis import Universe
from .trajectory import TrajectoryAnalyzer

def main():
    parser = argparse.ArgumentParser(
        description="Analyze metabolite-protein binding events from an MD trajectory."
    )
    parser.add_argument("topology", help="Topology file (e.g., PDB, PSF, GRO)")
    parser.add_argument("trajectory", help="Trajectory file (e.g., DCD, XTC, TRR)")
    parser.add_argument("--metabolite-selection", default=None,
                        help="MDAnalysis atom selection string for metabolites")
    parser.add_argument("--protein-selection", default=None,
                        help="MDAnalysis atom selection string for proteins")
    parser.add_argument("--cutoff", type=float, default=4.0,
                        help="Distance cutoff in Angstrom")
    parser.add_argument("--start", type=int, default=0,
                        help="Start frame")
    parser.add_argument("--stop", type=int, default=None,
                        help="Stop frame")
    parser.add_argument("--step", type=int, default=1,
                        help="Frame step size")
    parser.add_argument("--output", default="residues.pkl",
                        help="Output pickle file for molecule-level results")

    args = parser.parse_args()

    if not args.metabolite_selection:
        metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        metabolite_selection = 'resname ATP'
    if not args.protein_selection:
        protein_selection = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
    else:
        protein_selection = 'protein'

    # Load universe
    u = Universe(args.topology, args.trajectory)
    metabolites = u.select_atoms(metabolite_selection)
    proteins = u.select_atoms(protein_selection)

    # Run trajectory analysis
    analyzer = TrajectoryAnalyzer(u, metabolites, proteins,
                                  cutoff=args.cutoff,
                                  start=args.start,
                                  stop=args.stop,
                                  step=args.step)
    residues = analyzer.analyze()

    # Save molecule-level output
    with open(args.output, "wb") as f:
        pickle.dump(residues, f)

    print(f"Analysis complete. Results saved to {args.output}")

# Allow command-line execution via python -m binding_analysis
if __name__ == "__main__":
    main()
