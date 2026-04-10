import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm


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

def _worker_protein_chunk(args):
    """
    Worker function for processing a chunk of frames.
    Rebuilds the Universe inside the worker process.
    """
    topology, trajectory, frame_indices, prot_sel, r_max = args

    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    atoms = u.select_atoms(prot_sel)
    nodes = build_nodes(atoms)
    atom_molnums = atoms.molnums

    results = []
    for frame_idx in frame_indices:
        u.trajectory[frame_idx]
        frame_data = analyse_frame(atoms, atom_molnums, nodes, r_max)
        frame_data['frame'] = u.trajectory.ts.frame
        frame_data['time'] = u.trajectory.ts.time
        results.append(frame_data)
    return results


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
    import MDAnalysis as mda
    from concurrent.futures import ProcessPoolExecutor

    u = mda.Universe(topology, trajectory)
    stop = stop if stop is not None else len(u.trajectory)

    frame_indices = list(range(start, stop, step))
    chunks = [frame_indices[i:i + chunk_size]
              for i in range(0, len(frame_indices), chunk_size)]

    worker_args = [
        (topology, trajectory, chunk, prot_sel, r_max)
        for chunk in chunks
    ]

    all_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for chunk_results in tqdm(executor.map(_worker_protein_chunk, worker_args),
                                  total=len(chunks), desc="Analysing trajectory (parallel)"):
            all_results.extend(chunk_results)

    all_results.sort(key=lambda x: x['frame'])
    return all_results
