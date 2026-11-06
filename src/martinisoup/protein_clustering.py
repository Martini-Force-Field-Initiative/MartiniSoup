import numpy as np
from scipy.spatial import cKDTree
from MDAnalysis.analysis.base import ResultsGroup

from .analysis_base import ParallelAnalysisBase, n_parts_for_chunk_size


def build_nodes(atoms):
    """
    Build a node list from a protein atom group.

    Each node represents one protein molecule and records its molecule ID
    and moltype. The ID is 1-based to match NetworkX node-link convention.

    Parameters
    ----------
    atoms : MDAnalysis.AtomGroup

    Returns
    -------
    list of dict
        [{'id': int, 'type': str}, ...]
    """
    atom_molnums = atoms.molnums
    return [
        {'id': int(m) + 1, 'type': atoms.moltypes[atoms.molnums == m][0]}
        for m in np.unique(atom_molnums)
    ]


def analyse_frame(atoms, atom_molnums, nodes, r_max):
    """
    Build a contact network for a single trajectory frame.

    Contacts are computed at atom level using a KDTree, then lifted to
    molecule level so that each unique protein-protein pair appears once.
    The returned dict is in NetworkX node-link format.

    Parameters
    ----------
    atoms : MDAnalysis.AtomGroup
        Protein atoms with positions set to the current frame.
    atom_molnums : np.ndarray
        Per-atom molecule indices (precomputed from topology).
    nodes : list of dict
        Node list as returned by build_nodes.
    r_max : float
        Contact distance cutoff in Angstrom.

    Returns
    -------
    dict
        NetworkX node-link graph dict (without 'frame' / 'time' keys).
    """
    atoms.wrap()
    tree = cKDTree(atoms.positions, boxsize=atoms.dimensions[:3])
    pairs = tree.query_pairs(r_max)

    frame_contacts = {
        tuple(sorted((int(atom_molnums[i]), int(atom_molnums[j]))))
        for i, j in pairs
        if atom_molnums[i] != atom_molnums[j]
    }

    links = [{'source': src, 'target': tgt} for src, tgt in frame_contacts]

    return {
        'directed':   False,
        'multigraph': False,
        'graph':      {},
        'nodes':      nodes,
        'links':      links,
    }


class ProteinContactNetwork(ParallelAnalysisBase):
    """
    Per-frame protein-protein contact networks.

    For every frame, protein atoms within ``r_max`` of each other are found
    with a periodic KDTree and lifted to molecule level, giving one NetworkX
    node-link graph per frame (see :func:`analyse_frame`). In a parallel run
    each part builds the graphs for its own frames and the parts' lists are
    concatenated in frame order.

    Parameters
    ----------
    atoms : MDAnalysis.AtomGroup
        Protein atoms. Requires ``molnums`` (and ``moltypes`` if ``nodes`` is
        not given), e.g. from a GROMACS .tpr.
    r_max : float
        Contact distance cutoff in Angstrom (default: 6.0).
    nodes : list of dict, optional
        Node list as returned by :func:`build_nodes`. Built from ``atoms``
        if not given.
    **kwargs
        Passed to AnalysisBase (e.g. verbose).

    Results
    -------
    results.frames : list of dict
        Node-link graph dicts with added 'frame' and 'time' keys, in frame
        order.
    results.nodes : list of dict
        The node list shared by every frame.

    Examples
    --------
    >>> net = ProteinContactNetwork(u.select_atoms(prot_sel), r_max=6.0)
    >>> net.run(backend='multiprocessing', n_workers=4)
    >>> net.results.frames
    """

    def __init__(self, atoms, r_max=6.0, nodes=None, **kwargs):
        super().__init__(atoms.universe.trajectory, **kwargs)
        self._atoms = atoms
        self._atom_molnums = atoms.molnums
        self._nodes = build_nodes(atoms) if nodes is None else nodes
        self._r_max = r_max

    def _prepare(self):
        self.results.frames = []

    def _single_frame(self):
        frame_data = analyse_frame(
            self._atoms, self._atom_molnums, self._nodes, self._r_max,
        )
        frame_data['frame'] = self._ts.frame
        frame_data['time'] = self._ts.time
        self.results.frames.append(frame_data)

    def _get_aggregator(self):
        return ResultsGroup(lookup={'frames': ResultsGroup.flatten_sequence})

    def _conclude(self):
        self.results.nodes = self._nodes


def analyse_trajectory(u, atoms, nodes, r_max, start=0, stop=None, step=1):
    """
    Run per-frame protein contact network analysis over a trajectory slice
    (serial).

    Thin wrapper around :class:`ProteinContactNetwork` run with the serial
    backend.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Unused beyond consistency with other analyses; the trajectory is taken
        from ``atoms.universe``.
    atoms : MDAnalysis.AtomGroup
        Protein atoms.
    nodes : list of dict
        Node list as returned by build_nodes.
    r_max : float
        Contact distance cutoff in Angstrom.
    start, stop, step : int
        Frame range, as in ``u.trajectory[start:stop:step]``.

    Returns
    -------
    list of dict
        Each entry is a node-link graph dict with added 'frame' and 'time'
        keys, in frame order.
    """
    net = ProteinContactNetwork(atoms, r_max=r_max, nodes=nodes, verbose=True)
    net.run(start=start, stop=stop, step=step)
    return net.results.frames


def analyse_trajectory_parallel(
    topology,
    trajectory,
    prot_sel,
    r_max=6.0,
    start=0,
    stop=None,
    step=1,
    n_workers=4,
    chunk_size=100,
):
    """
    Parallel per-frame protein contact network analysis.

    Thin wrapper around :class:`ProteinContactNetwork` run with MDAnalysis's
    ``multiprocessing`` backend. The frame range is split into contiguous
    parts of roughly ``chunk_size`` frames, distributed over ``n_workers``
    processes.

    Parameters
    ----------
    topology, trajectory : str
        Paths to topology and trajectory files.
    prot_sel : str
        MDAnalysis selection string for proteins.
    r_max : float
        Contact distance cutoff in Angstrom.
    start, stop, step : int
        Frame range, as in ``u.trajectory[start:stop:step]``.
    n_workers : int
        Number of worker processes.
    chunk_size : int
        Approximate frames per part.

    Returns
    -------
    list of dict
        Per-frame node-link graph dicts, in frame order.
    """
    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    n_parts = n_parts_for_chunk_size(u.trajectory, start, stop, step, chunk_size)

    net = ProteinContactNetwork(u.select_atoms(prot_sel), r_max=r_max)
    net.run(
        start=start, stop=stop, step=step,
        backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
    )
    return net.results.frames
