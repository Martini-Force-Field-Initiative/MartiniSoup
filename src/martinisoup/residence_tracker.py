from tqdm import tqdm
from scipy.spatial import cKDTree as KDTree
from .data_structures import BindingState, ResidenceRegistry


# --------------------------------------------------------------------------- #
# Parallel helpers
# --------------------------------------------------------------------------- #

def _compute_bound_atoms_chunk(args):
    """
    Worker: jump to each frame in *frame_indices*, run the KDTree query, and
    return a list of (frame_idx, stamp, frozenset_of_bound_metabolite_atom_indices).

    Rebuilds the Universe inside the worker process so the main process does
    not need to share any MDAnalysis state.
    """
    (topology, trajectory, prot_sel, metab_sel, cutoff, frame_indices, use_time) = args

    from MDAnalysis import Universe, transformations

    u = Universe(topology, trajectory)
    ag = u.atoms
    u.trajectory.add_transformations(transformations.wrap(ag))

    proteins = u.select_atoms(prot_sel)
    metabolites = u.select_atoms(metab_sel)

    results = []
    for frame_idx in frame_indices:
        u.trajectory[frame_idx]
        stamp = u.trajectory.ts.time if use_time else u.trajectory.ts.frame

        protein_tree  = KDTree(proteins.positions,   boxsize=u.dimensions[:3])
        metabolite_tree = KDTree(metabolites.positions, boxsize=u.dimensions[:3])
        sdm = protein_tree.sparse_distance_matrix(metabolite_tree, cutoff)

        bound_met_atoms = frozenset(pair[1] for pair in sdm.keys())
        results.append((frame_idx, stamp, bound_met_atoms))

    return results


def track_parallel(
    topology,
    trajectory,
    metab_sel,
    prot_sel,
    cutoff=4.0,
    start=0,
    stop=None,
    step=1,
    n_workers=4,
    chunk_size=100,
    use_time=True,
):
    """
    Two-phase parallel residence-time analysis.

    Phase 1 (parallel)
        Workers compute, for each frame, the set of metabolite atom indices
        within *cutoff* of any protein atom.  Each frame is independent so
        this is embarrassingly parallel.

    Phase 2 (sequential, cheap)
        The bound-sets are reassembled in frame order and fed into the same
        BindingState machine used by BindingEventTracker.  Because Phase 2
        sees the full ordered sequence, binding events that span chunk
        boundaries are recorded correctly.

    Returns
    -------
    dict
        {moltype_name: [durations]} — identical structure to
        BindingEventTracker.track().
    """
    from MDAnalysis import Universe, transformations
    from concurrent.futures import ProcessPoolExecutor

    # ------------------------------------------------------------------ #
    # Setup: moltype mapping and atom metadata (done once in main process)
    # ------------------------------------------------------------------ #
    u = Universe(topology, trajectory)
    ag = u.atoms
    u.trajectory.add_transformations(transformations.wrap(ag))

    n_frames = len(u.trajectory)
    stop = stop if stop is not None else n_frames

    metabolites = u.select_atoms(metab_sel)

    unique = {}
    moltype_table = []
    moltype_ids = []
    for mt in metabolites.moltypes:
        if mt not in unique:
            unique[mt] = len(moltype_table)
            moltype_table.append(mt)
        moltype_ids.append(unique[mt])

    molnums  = list(metabolites.molnums)
    n_atoms  = len(metabolites)

    # ------------------------------------------------------------------ #
    # Phase 1: parallel KDTree queries
    # ------------------------------------------------------------------ #
    frame_indices = list(range(start, stop, step))
    chunks = [frame_indices[i:i + chunk_size]
              for i in range(0, len(frame_indices), chunk_size)]

    worker_args = [
        (topology, trajectory, prot_sel, metab_sel, cutoff, chunk, use_time)
        for chunk in chunks
    ]

    all_frame_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for result in tqdm(executor.map(_compute_bound_atoms_chunk, worker_args),
                           total=len(chunks), desc="Computing contacts (parallel)"):
            all_frame_results.extend(result)

    # Sort so Phase 2 sees frames in trajectory order
    all_frame_results.sort(key=lambda x: x[0])

    # ------------------------------------------------------------------ #
    # Phase 2: sequential state machine (cheap — no geometry)
    # ------------------------------------------------------------------ #
    tracker = {i: BindingState(moltype_id=moltype_ids[i]) for i in range(n_atoms)}

    last_stamp = None
    for _, stamp, bound_met_atoms in tqdm(all_frame_results,
                                                   desc="Running state machine"):
        last_stamp = stamp
        for atom_idx, state in tracker.items():
            state.update_bound(atom_idx in bound_met_atoms, stamp)

    final_stamp = last_stamp if last_stamp is not None else 0
    for state in tracker.values():
        state.finalize(final_stamp)

    return ResidenceRegistry(
        tracker=tracker,
        molnums=molnums,
        moltype_table={i: name for i, name in enumerate(moltype_table)},
        moltype_ids=list(moltype_ids),
    ).get_durations_by_type()

class BindingEventTracker:
    def __init__(self, u, metabolites, proteins, cutoff=4.0,
                 start=0, stop=None, step=1, use_time=True):
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

    def track(self):
        n_atoms = len(self.metabolites)
        # initialize tracker with moltype ids
        tracker = {i: BindingState(moltype_id=self.moltype_ids[i]) for i in range(n_atoms)}

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
        residues = ResidenceRegistry(
            tracker=tracker,
            molnums=list(self.metabolites.molnums),   # ensure serializable list-like
            moltype_table={i: name for i, name in enumerate(self.moltype_table)},
            moltype_ids=list(self.moltype_ids),
        ).get_durations_by_type()

        return residues
