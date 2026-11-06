from tqdm import tqdm
from scipy.spatial import cKDTree as KDTree
from MDAnalysis import Universe
from MDAnalysis.analysis.base import ResultsGroup
from MDAnalysis.exceptions import NoDataError

from .analysis_base import ParallelAnalysisBase, n_parts_for_chunk_size
from .data_structures import BindingState, ResidenceRegistry


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
# Per-frame contact query
# --------------------------------------------------------------------------- #

def _bound_metabolite_atoms(proteins, metabolites, box, cutoff):
    """
    Indices (into ``metabolites``) of metabolite atoms within ``cutoff`` of
    any protein atom, using periodic KDTrees.

    Positions are wrapped into the unit cell first (cKDTree's ``boxsize``
    requires it). This uses the same calculation as the ``wrap``
    transformation, but leaves the Universe untouched.
    """
    boxsize = box[:3]
    protein_tree = KDTree(proteins.wrap(compound='atoms', inplace=False), boxsize=boxsize)
    metabolite_tree = KDTree(metabolites.wrap(compound='atoms', inplace=False), boxsize=boxsize)
    sdm = protein_tree.sparse_distance_matrix(metabolite_tree, cutoff)
    return frozenset(pair[1] for pair in sdm.keys())


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

    The fold is inherently sequential and stateful (a binding run may span
    many frames, and therefore many parallel parts), so it is never
    parallelised: ResidenceTracker runs it once in ``_conclude``, in the main
    process, after the per-frame contact sets from all parts have been
    concatenated in frame order.
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


class ResidenceTracker(ParallelAnalysisBase):
    """
    Metabolite residence times on proteins.

    Phase 1 (per frame, parallelisable)
        ``_single_frame`` finds the metabolite atoms within ``cutoff`` of any
        protein atom. Each frame is independent, so parts of the trajectory
        can be processed by separate workers; their per-frame contact sets
        are concatenated in frame order.

    Phase 2 (sequential, in ``_conclude``)
        The ordered contact sets are folded through a per-atom
        bound/unbound state machine to give binding durations. Because this
        sees the full ordered sequence, binding events that span part
        boundaries are recorded as single durations.

    Parameters
    ----------
    metabolites : MDAnalysis.AtomGroup
        Metabolite atoms. Uses ``moltypes``/``molnums`` if present (e.g. from
        a .tpr), falling back to ``resnames``/``resindices``.
    proteins : MDAnalysis.AtomGroup
    cutoff : float
        Contact distance cutoff in Angstrom (default: 4.0).
    use_time : bool
        Stamp events with ``ts.time`` (True, default) or ``ts.frame``
        (False); durations are differences of these stamps.
    **kwargs
        Passed to AnalysisBase (e.g. verbose).

    Results
    -------
    results.durations : dict
        {moltype_name: [durations]}

    Examples
    --------
    >>> rt = ResidenceTracker(metabolites, proteins, cutoff=5.0)
    >>> rt.run(backend='multiprocessing', n_workers=4)
    >>> rt.results.durations
    """

    def __init__(self, metabolites, proteins, cutoff=4.0, use_time=True, **kwargs):
        super().__init__(metabolites.universe.trajectory, **kwargs)
        self._metabolites = metabolites
        self._proteins = proteins
        self._cutoff = cutoff
        self._use_time = use_time

        self._moltype_table, self._moltype_ids = _build_moltype_mapping(metabolites)
        self._molnums = list(_molnums(metabolites))

    def _prepare(self):
        self.results.frame_contacts = []

    def _single_frame(self):
        ts = self._ts
        stamp = ts.time if self._use_time else ts.frame
        bound = _bound_metabolite_atoms(
            self._proteins, self._metabolites, ts.dimensions, self._cutoff,
        )
        self.results.frame_contacts.append((ts.frame, stamp, bound))

    def _get_aggregator(self):
        return ResultsGroup(lookup={'frame_contacts': ResultsGroup.flatten_sequence})

    def _conclude(self):
        # Parts are already concatenated in frame order, but the state
        # machine's correctness depends on strict ordering, so make that
        # guarantee explicit rather than implicit.
        frame_contacts = sorted(self.results.frame_contacts, key=lambda x: x[0])
        self.results.durations = _run_state_machine(
            frame_contacts, self._moltype_ids, len(self._metabolites),
            self._molnums, self._moltype_table,
            desc="Running state machine", total=len(frame_contacts),
        )
        # The per-frame contact sets are only needed for the fold.
        del self.results['frame_contacts']


# --------------------------------------------------------------------------- #
# Function wrappers (used by the residence-times command)
# --------------------------------------------------------------------------- #

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
    Single-process residence-time analysis.

    Thin wrapper around :class:`ResidenceTracker` run with the serial
    backend.

    Returns
    -------
    dict
        {moltype_name: [durations]}
    """
    u = Universe(topology, trajectory)
    rt = ResidenceTracker(
        u.select_atoms(metab_sel), u.select_atoms(prot_sel),
        cutoff=cutoff, use_time=use_time, verbose=True,
    )
    rt.run(start=start, stop=stop, step=step)
    return rt.results.durations


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
    Parallel residence-time analysis.

    Thin wrapper around :class:`ResidenceTracker` run with MDAnalysis's
    ``multiprocessing`` backend. The per-frame contact query is split into
    contiguous parts of roughly ``chunk_size`` frames over ``n_workers``
    processes; the state machine then runs once, in order, in the main
    process, so results are identical to :func:`track_serial`.

    Returns
    -------
    dict
        {moltype_name: [durations]}
    """
    u = Universe(topology, trajectory)
    n_parts = n_parts_for_chunk_size(u.trajectory, start, stop, step, chunk_size)

    rt = ResidenceTracker(
        u.select_atoms(metab_sel), u.select_atoms(prot_sel),
        cutoff=cutoff, use_time=use_time,
    )
    rt.run(
        start=start, stop=stop, step=step,
        backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
    )
    return rt.results.durations
