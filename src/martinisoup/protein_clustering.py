import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

from .parallel import map_trajectory_parallel


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


def analyse_trajectory(u, atoms, nodes, r_max, start=0, stop=None, step=1):
    """
    Run per-frame protein contact network analysis over a trajectory slice.

    Parameters
    ----------
    u : MDAnalysis.Universe
    atoms : MDAnalysis.AtomGroup
        Protein atoms.
    nodes : list of dict
        Node list as returned by build_nodes.
    r_max : float
        Contact distance cutoff in Angstrom.
    start : int
        First frame index.
    stop : int or None
        Last frame index (exclusive). None means end of trajectory.
    step : int
        Frame stride.

    Returns
    -------
    list of dict
        Each entry is a node-link graph dict with added 'frame' and 'time' keys.
    """
    atom_molnums = atoms.molnums
    results = []
    for ts in tqdm(u.trajectory[start:stop:step]):
        frame_data = analyse_frame(atoms, atom_molnums, nodes, r_max)
        frame_data['frame'] = ts.frame
        frame_data['time'] = ts.time
        results.append(frame_data)
    return results


# --------------------------------------------------------------------------- #
# Parallel version (multiprocessing)
# --------------------------------------------------------------------------- #

def _protein_setup(u, prot_sel, r_max):
    """Build the per-worker context: protein atoms, molecule numbers, nodes."""
    atoms = u.select_atoms(prot_sel)
    nodes = build_nodes(atoms)
    return atoms, atoms.molnums, nodes, r_max


def _protein_per_frame(u, context, ts):
    """Build one frame's contact network."""
    atoms, atom_molnums, nodes, r_max = context
    frame_data = analyse_frame(atoms, atom_molnums, nodes, r_max)
    frame_data['frame'] = ts.frame
    frame_data['time'] = ts.time
    return frame_data


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

    Each frame is independent, so the trajectory slice is split into chunks
    and dispatched to worker processes. Results are collected and returned
    sorted by frame index.

    Parameters
    ----------
    topology, trajectory : str
        Paths to topology and trajectory files.
    prot_sel : str
        MDAnalysis selection string for proteins.
    r_max : float
        Contact distance cutoff in Angstrom.
    start : int
        First frame index.
    stop : int or None
        Last frame index (exclusive). None means end of trajectory.
    step : int
        Frame stride.
    n_workers : int
        Number of worker processes.
    chunk_size : int
        Frames per chunk.

    Returns
    -------
    list of dict
        Per-frame node-link graph dicts, sorted by frame index.
    """
    all_results = map_trajectory_parallel(
        topology, trajectory, _protein_setup, _protein_per_frame,
        setup_args=(prot_sel, r_max),
        start=start, stop=stop, step=step,
        n_workers=n_workers, chunk_size=chunk_size,
        desc="Analysing trajectory (parallel)",
    )
    all_results.sort(key=lambda x: x['frame'])
    return all_results
