from tqdm import tqdm
from scipy.spatial import cKDTree as KDTree
from .data_structures import BindingState, ResidenceRegistry
from .parallel import map_trajectory_parallel
from MDAnalysis import Universe, transformations
from MDAnalysis.exceptions import NoDataError


def _moltypes(ag):
    """Return per-atom moltype labels, falling back to resnames if unavailable."""
    try:
        return ag.moltypes
    except NoDataError:
        return ag.resnames


def _molnums(ag):
    """Return per-atom molecule indices, falling back to resindices if unavailable."""
    try:
        return ag.molnums
    except NoDataError:
        return ag.resindices


def _build_moltype_mapping(metabolites):
    """
    Build a stable moltype -> integer id mapping and the per-atom id list.

    Returns
    -------
    moltype_table : dict
        {moltype_id: moltype_name}
    moltype_ids : list of int
        Per-atom moltype id, in atom order.
    """
    unique = {}
    names = []
    moltype_ids = []
    for mt in _moltypes(metabolites):
        if mt not in unique:
            unique[mt] = len(names)
            names.append(mt)
        moltype_ids.append(unique[mt])
    return {i: name for i, name in enumerate(names)}, moltype_ids


# --------------------------------------------------------------------------- #
# Per-frame contact query (shared by serial and parallel paths)
# --------------------------------------------------------------------------- #

def _residence_setup(u, prot_sel, metab_sel, cutoff, use_time):
    """Build the per-worker context: wrapped atoms, selections, cutoff, stamp mode."""
    ag = u.atoms
    u.trajectory.add_transformations(transformations.wrap(ag))
    proteins = u.select_atoms(prot_sel)
    metabolites = u.select_atoms(metab_sel)
    return proteins, metabolites, cutoff, use_time


def _residence_per_frame(u, context, ts):
    """
    Find metabolite atoms within *cutoff* of any protein atom for one frame.

    Returns (frame_idx, stamp, frozenset_of_bound_metabolite_atom_indices).
    """
    proteins, metabolites, cutoff, use_time = context
    stamp = ts.time if use_time else ts.frame

    protein_tree = KDTree(proteins.positions, boxsize=u.dimensions[:3])
    metabolite_tree = KDTree(metabolites.positions, boxsize=u.dimensions[:3])
    sdm = protein_tree.sparse_distance_matrix(metabolite_tree, cutoff)

    bound_met_atoms = frozenset(pair[1] for pair in sdm.keys())
    return ts.frame, stamp, bound_met_atoms


# --------------------------------------------------------------------------- #
# State machine fold (shared by serial and parallel paths)
# --------------------------------------------------------------------------- #

def _run_state_machine(
    frame_results, moltype_ids, n_atoms, molnums, moltype_table,
    desc="Running state machine", total=None,
):
    """
    Fold an ordered sequence of (frame_idx, stamp, bound_atom_indices) triples
    into final per-moltype residence durations.

    This is the same computation whether the (frame_idx, stamp, bound_atoms)
    triples were produced by stepping through the trajectory directly
    (track_serial) or by a parallel geometry pass followed by reassembly in
    frame order (track_parallel) — the fold itself is inherently sequential
    and stateful (a binding run may span many frames), so it is never
    parallelised, only fed from either source.
    """
    tracker = {i: BindingState(moltype_id=moltype_ids[i]) for i in range(n_atoms)}

    last_stamp = None
    for _, stamp, bound_met_atoms in tqdm(frame_results, desc=desc, total=total):
        last_stamp = stamp
        for atom_idx, state in tracker.items():
            state.update_bound(atom_idx in bound_met_atoms, stamp)

    final_stamp = last_stamp if last_stamp is not None else 0
    for state in tracker.values():
        state.finalize(final_stamp)

    return ResidenceRegistry(
        tracker=tracker,
        molnums=molnums,
        moltype_table=moltype_table,
        moltype_ids=moltype_ids,
    ).get_durations_by_type()


def track_serial(
    topology,
    trajectory,
    metab_sel,
    prot_sel,
    cutoff=4.0,
    start=0,
    stop=None,
    step=1,
    use_time=True,
):
    """
    Single-threaded residence-time analysis.

    Returns
    -------
    dict
        {moltype_name: [durations]}
    """
    u = Universe(topology, trajectory)
    context = _residence_setup(u, prot_sel, metab_sel, cutoff, use_time)
    metabolites = context[1]

    stop = stop if stop is not None else len(u.trajectory)

    moltype_table, moltype_ids = _build_moltype_mapping(metabolites)
    n_atoms = len(metabolites)
    molnums = list(_molnums(metabolites))

    frame_indices = range(start, stop, step)
    frame_results = (
        _residence_per_frame(u, context, u.trajectory[i]) for i in frame_indices
    )

    return _run_state_machine(
        frame_results, moltype_ids, n_atoms, molnums, moltype_table,
        desc="Analyzing trajectory", total=len(frame_indices),
    )


# --------------------------------------------------------------------------- #
# Parallel version (multiprocessing)
# --------------------------------------------------------------------------- #

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
        this is embarrassingly parallel, and uses the same per-frame query
        (_residence_per_frame) as track_serial.

    Phase 2 (sequential, cheap)
        The bound-sets are reassembled in frame order and folded through
        _run_state_machine — the exact same state-machine fold track_serial
        uses, just fed from a pre-computed list instead of a live trajectory
        iterator. Because Phase 2 sees the full ordered sequence, binding
        events that span chunk boundaries are recorded correctly.

    Returns
    -------
    dict
        {moltype_name: [durations]}
    """
    # ------------------------------------------------------------------ #
    # Setup: moltype mapping and atom metadata (done once in main process)
    # ------------------------------------------------------------------ #
    u = Universe(topology, trajectory)
    metabolites = u.select_atoms(metab_sel)

    moltype_table, moltype_ids = _build_moltype_mapping(metabolites)
    n_atoms = len(metabolites)
    molnums = list(_molnums(metabolites))

    # ------------------------------------------------------------------ #
    # Phase 1: parallel KDTree queries
    # ------------------------------------------------------------------ #
    frame_results = map_trajectory_parallel(
        topology, trajectory, _residence_setup, _residence_per_frame,
        setup_args=(prot_sel, metab_sel, cutoff, use_time),
        start=start, stop=stop, step=step,
        n_workers=n_workers, chunk_size=chunk_size,
        desc="Computing contacts (parallel)",
    )

    # Defensive: map_trajectory_parallel already returns results in frame
    # order, but Phase 2 correctness depends on strict ordering, so make
    # that guarantee explicit rather than implicit.
    frame_results.sort(key=lambda x: x[0])

    # ------------------------------------------------------------------ #
    # Phase 2: sequential state machine (cheap — no geometry)
    # ------------------------------------------------------------------ #
    return _run_state_machine(
        frame_results, moltype_ids, n_atoms, molnums, moltype_table,
        desc="Running state machine", total=len(frame_results),
    )
