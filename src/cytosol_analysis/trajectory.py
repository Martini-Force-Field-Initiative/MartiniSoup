from tqdm import tqdm
from scipy.spatial import cKDTree as KDTree
from .data_structures import BoundState, MetaboliteResidences

class TrajectoryAnalyzer:
    def __init__(self, u, metabolites, proteins, cutoff=4.0,
                 start=0, stop=None, step=1, use_time=False):
        """
        use_time: if True, durations are in physical time (ts.time),
                  otherwise durations are in frame counts (ts.frame).
        """
        self.u = u
        self.metabolites = metabolites
        self.proteins = proteins
        self.cutoff = cutoff
        self.start = start
        self.stop = stop if stop is not None else len(u.trajectory)
        self.step = step
        self.use_time = use_time

        # build moltype mapping (per-atom)
        unique = {}
        moltype_table = []
        moltype_ids = []
        for mt in metabolites.moltypes:          # MDAnalysis per-atom property
            if mt not in unique:
                unique[mt] = len(moltype_table)
                moltype_table.append(mt)
            moltype_ids.append(unique[mt])

        self.moltype_ids = moltype_ids     # list, len = n_metabolite_atoms
        self.moltype_table = moltype_table # list: id -> name

    def run(self):
        n_atoms = len(self.metabolites)
        # initialize tracker with moltype ids
        tracker = {i: BoundState(moltype_id=self.moltype_ids[i]) for i in range(n_atoms)}

        # iterate trajectory
        last_stamp = None
        for ts in tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Analyzing trajectory"):
            # choose stamp to record (frame index or time)
            stamp = ts.time if self.use_time else ts.frame
            last_stamp = stamp

            # build periodic KDTree for protein and metabolite positions
            # note: KDTree(..., boxsize=...) expects boxsize as scalar or array -> use u.dimensions[:3]
            protein_tree = KDTree(self.proteins.positions, boxsize=self.u.dimensions[:3])
            metabolite_tree = KDTree(self.metabolites.positions, boxsize=self.u.dimensions[:3])

            # sparse distance matrix: keys are (i_from_protein, j_from_metabolite)
            sdm = protein_tree.sparse_distance_matrix(metabolite_tree, self.cutoff)

            # collect metabolite atom indices that are within cutoff of any protein atom
            # sdm is a dict-like: keys -> (i_protein, j_metabolite)
            bound_met_atoms = {pair[1] for pair in sdm.keys()}  # set of metabolite indices

            # update each atom's BoundState
            for atom_idx, state in tracker.items():
                is_bound = (atom_idx in bound_met_atoms)
                state.update_bound(is_bound, stamp)

        # finalize any open events using the last stamp seen
        if last_stamp is None:
            # no frames processed; nothing to finalize
            final_stamp = 0
        else:
            final_stamp = last_stamp

        for state in tracker.values():
            state.finalize(final_stamp)

        # Build MetaboliteResidences, including moltype_table and moltype_ids
        residues = MetaboliteResidences(
            tracker=tracker,
            molnums=list(self.metabolites.molnums),   # ensure serializable list-like
            moltype_table={i: name for i, name in enumerate(self.moltype_table)},
            moltype_ids=list(self.moltype_ids),
        ).make_type_agg_named()

        return residues
