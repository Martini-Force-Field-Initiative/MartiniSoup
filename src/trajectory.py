from tqdm import tqdm
from scipy.spatial import cKDTree as KDTree
from data_structures import BoundState, MetaboliteResidences
import pickle

class TrajectoryAnalyzer:
    """
    Analyze MD trajectories to detect metabolite-protein binding events.

    Parameters
    ----------
    u : MDAnalysis.Universe
        The MDAnalysis universe containing trajectory and topology.
    metabolites : MDAnalysis AtomGroup
        AtomGroup for metabolites.
    proteins : MDAnalysis AtomGroup
        AtomGroup for proteins.
    cutoff : float
        Distance cutoff (Angstrom) to define binding.
    start : int
        First frame to analyze.
    stop : int
        Last frame to analyze.
    step : int
        Frame step size.
    """

    def __init__(self, u, metabolites, proteins, cutoff=4.0, start=None, stop=None, step=1):
        self.u = u
        self.metabolites = metabolites
        self.proteins = proteins
        self.cutoff = cutoff
        self.start = start if start is not None else int(len(u.trajectory) / 2)
        self.stop = stop if stop is not None else len(u.trajectory)
        self.step = step

        self.n_metabolites = len(metabolites)
        # Assuming each metabolite atom has a moltype_id and molecule index
        self.moltypes = [getattr(atom, 'moltype_id', None) for atom in metabolites]
        self.molnums = [getattr(atom, 'molnum', None) for atom in metabolites]

    def analyze(self, savepath='.'):
        """
        Run the trajectory analysis and return atom-level BoundState tracker
        and a MetaboliteResidences object.
        """
        # Initialize tracker
        tracker = {i: BoundState(moltype_id=mt) for i, mt in enumerate(self.moltypes)}

        for ts in tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Analyzing trajectory"):
            # Build periodic KDTree for protein and metabolite positions
            protein_tree = KDTree(self.proteins.positions, boxsize=self.u.dimensions[:3])
            metabolite_tree = KDTree(self.metabolites.positions, boxsize=self.u.dimensions[:3])

            # Compute sparse distance matrix with cutoff
            sdm = protein_tree.sparse_distance_matrix(metabolite_tree, self.cutoff)

            # Convert sparse distance matrix keys to metabolite atom indices
            bound_atoms = set(j for (_, j), dist in sdm.items())

            # Update tracker
            for atom_idx, state in tracker.items():
                if atom_idx in bound_atoms:
                    if not state.bound:
                        # New binding event
                        state.bound = True
                        state.start = ts.frame
                else:
                    if state.bound:
                        # Event ended
                        duration = ts.frame - state.start
                        state.durations.append(duration)
                        state.bound = False
                        state.start = None

        # Ensure any remaining bound states are closed at the end
        for state in tracker.values():
            if state.bound:
                duration = ts.frame - state.start
                state.durations.append(duration)
                state.bound = False
                state.start = None

        # Construct MetaboliteResidences object
        residues = MetaboliteResidences(tracker, self.molnums, moltype_table=None)
        if savepath:
            with open(savepath, 'wb') as f:
                pickle.dump(residues, f)
        return residues
