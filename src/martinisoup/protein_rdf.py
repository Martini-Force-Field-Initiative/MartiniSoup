import numpy as np
from tqdm import tqdm
from MDAnalysis.lib.distances import self_distance_array


def compute_rdf(u, atoms, r_max=100.0, n_bins=200, start=0, stop=None, step=1):
    """
    Compute the protein centre-of-mass RDF over a trajectory slice.

    Each protein molecule's centre of mass is computed per frame, then all
    pairwise distances are histogrammed. The RDF is normalised by the mean
    box volume and the number of protein pairs.

    Parameters
    ----------
    u : MDAnalysis.Universe
    atoms : MDAnalysis.AtomGroup
        Protein atoms. Molecules are inferred via atoms.split('molecule').
    r_max : float
        Upper distance limit in Angstrom (default: 100.0).
    n_bins : int
        Number of histogram bins (default: 200).
    start : int
        First frame index (default: 0).
    stop : int or None
        Last frame index, exclusive. None means end of trajectory.
    step : int
        Frame stride (default: 1).

    Returns
    -------
    dict with keys:
        'r'          : np.ndarray, bin midpoints (Angstrom)
        'gr'         : np.ndarray, g(r) values
        'n_proteins' : int, number of protein molecules
        'mean_vol'   : float, mean box volume over analysed frames (Angstrom^3)
    """
    protein_molecules = atoms.split('molecule')
    n_proteins = len(protein_molecules)

    bins    = np.linspace(0, r_max, n_bins + 1)
    counts  = np.zeros(n_bins, dtype=np.float64)
    vol_acc = 0.0
    n_frames = 0

    for ts in tqdm(u.trajectory[start:stop:step], desc="Computing protein RDF"):
        coms    = np.array([mol.center_of_mass(wrap=True) for mol in protein_molecules])
        counts += np.histogram(self_distance_array(coms, box=ts.dimensions), bins=bins)[0]
        vol_acc += np.prod(ts.dimensions[:3])
        n_frames += 1

    mean_vol  = vol_acc / n_frames
    r_mid     = 0.5 * (bins[:-1] + bins[1:])
    shell_vol = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)
    gr        = (2 * mean_vol) / (n_proteins * (n_proteins - 1)) * counts / (shell_vol * n_frames)

    return {'r': r_mid, 'gr': gr, 'n_proteins': n_proteins, 'mean_vol': mean_vol}


# --------------------------------------------------------------------------- #
# Parallel version (multiprocessing)
# --------------------------------------------------------------------------- #

def _worker_rdf_chunk(args):
    """
    Worker function for processing a chunk of frames.

    Returns partial histogram counts and accumulated box volume for the
    assigned frames. These are summed in the main process before normalisation.
    Rebuilds the Universe inside the worker process to avoid pickling issues.
    """
    topology, trajectory, frame_indices, prot_sel, r_max, n_bins = args

    import MDAnalysis as mda
    from MDAnalysis.lib.distances import self_distance_array

    u = mda.Universe(topology, trajectory)
    protein_molecules = u.select_atoms(prot_sel).split('molecule')
    n_proteins = len(protein_molecules)

    bins    = np.linspace(0, r_max, n_bins + 1)
    counts  = np.zeros(n_bins, dtype=np.float64)
    vol_acc = 0.0

    for frame_idx in frame_indices:
        u.trajectory[frame_idx]
        ts   = u.trajectory.ts
        coms = np.array([mol.center_of_mass(wrap=True) for mol in protein_molecules])
        counts  += np.histogram(self_distance_array(coms, box=ts.dimensions), bins=bins)[0]
        vol_acc += np.prod(ts.dimensions[:3])

    return {'counts': counts, 'vol_acc': vol_acc, 'n_proteins': n_proteins}


def compute_rdf_parallel(
    topology,
    trajectory,
    prot_sel,
    r_max=100.0,
    n_bins=200,
    start=0,
    stop=None,
    step=1,
    n_workers=4,
    chunk_size=100,
):
    """
    Compute the protein centre-of-mass RDF in parallel.

    The trajectory slice is split into frame-index chunks and dispatched to
    worker processes. Each worker accumulates partial histogram counts and
    box volumes; results are merged and normalised in the main process.

    Parameters
    ----------
    topology, trajectory : str
        Paths to topology and trajectory files.
    prot_sel : str
        MDAnalysis selection string for protein atoms.
    r_max : float
        Upper distance limit in Angstrom (default: 100.0).
    n_bins : int
        Number of histogram bins (default: 200).
    start : int
        First frame index (default: 0).
    stop : int or None
        Last frame index, exclusive. None means end of trajectory.
    step : int
        Frame stride (default: 1).
    n_workers : int
        Number of worker processes (default: 4).
    chunk_size : int
        Frames per chunk (default: 100).

    Returns
    -------
    dict with keys:
        'r'          : np.ndarray, bin midpoints (Angstrom)
        'gr'         : np.ndarray, g(r) values
        'n_proteins' : int, number of protein molecules
        'mean_vol'   : float, mean box volume over analysed frames (Angstrom^3)
    """
    import MDAnalysis as mda
    from concurrent.futures import ProcessPoolExecutor

    u    = mda.Universe(topology, trajectory)
    stop = stop if stop is not None else len(u.trajectory)

    frame_indices = list(range(start, stop, step))
    n_frames      = len(frame_indices)
    chunks        = [frame_indices[i:i + chunk_size]
                     for i in range(0, n_frames, chunk_size)]

    worker_args = [
        (topology, trajectory, chunk, prot_sel, r_max, n_bins)
        for chunk in chunks
    ]

    bins         = np.linspace(0, r_max, n_bins + 1)
    total_counts = np.zeros(n_bins, dtype=np.float64)
    total_vol    = 0.0
    n_proteins   = None

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for result in tqdm(executor.map(_worker_rdf_chunk, worker_args),
                           total=len(chunks), desc="Computing protein RDF (parallel)"):
            total_counts += result['counts']
            total_vol    += result['vol_acc']
            if n_proteins is None:
                n_proteins = result['n_proteins']

    mean_vol  = total_vol / n_frames
    r_mid     = 0.5 * (bins[:-1] + bins[1:])
    shell_vol = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)
    gr        = (2 * mean_vol) / (n_proteins * (n_proteins - 1)) * total_counts / (shell_vol * n_frames)

    return {'r': r_mid, 'gr': gr, 'n_proteins': n_proteins, 'mean_vol': mean_vol}
