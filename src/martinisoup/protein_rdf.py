import numpy as np
from tqdm import tqdm
from MDAnalysis.lib.distances import self_distance_array

from .parallel import map_trajectory_parallel


def _rdf_setup(u, prot_sel, r_max, n_bins):
    """Build the per-worker context: protein molecules and histogram bins."""
    protein_molecules = u.select_atoms(prot_sel).split('molecule')
    bins = np.linspace(0, r_max, n_bins + 1)
    return protein_molecules, bins


def _rdf_per_frame(u, context, ts):
    """Histogram one frame's protein centre-of-mass pairwise distances."""
    protein_molecules, bins = context
    coms = np.array([mol.center_of_mass(wrap=True) for mol in protein_molecules])
    counts = np.histogram(self_distance_array(coms, box=ts.dimensions), bins=bins)[0]
    return counts, np.prod(ts.dimensions[:3])


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
    bins = np.linspace(0, r_max, n_bins + 1)
    context = (protein_molecules, bins)

    counts = np.zeros(n_bins, dtype=np.float64)
    vol_acc = 0.0
    n_frames = 0

    for ts in tqdm(u.trajectory[start:stop:step], desc="Computing protein RDF"):
        frame_counts, frame_vol = _rdf_per_frame(u, context, ts)
        counts += frame_counts
        vol_acc += frame_vol
        n_frames += 1

    mean_vol  = vol_acc / n_frames
    r_mid     = 0.5 * (bins[:-1] + bins[1:])
    shell_vol = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)
    gr        = (2 * mean_vol) / (n_proteins * (n_proteins - 1)) * counts / (shell_vol * n_frames)

    return {'r': r_mid, 'gr': gr, 'n_proteins': n_proteins, 'mean_vol': mean_vol}


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
    worker processes via `map_trajectory_parallel`. Each frame's histogram
    counts and box volume are merged and normalised in the main process.

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

    u = mda.Universe(topology, trajectory)
    n_proteins = len(u.select_atoms(prot_sel).split('molecule'))

    per_frame_results = map_trajectory_parallel(
        topology, trajectory, _rdf_setup, _rdf_per_frame,
        setup_args=(prot_sel, r_max, n_bins),
        start=start, stop=stop, step=step,
        n_workers=n_workers, chunk_size=chunk_size,
        desc="Computing protein RDF (parallel)",
    )

    bins = np.linspace(0, r_max, n_bins + 1)
    total_counts = np.zeros(n_bins, dtype=np.float64)
    total_vol = 0.0
    for frame_counts, frame_vol in per_frame_results:
        total_counts += frame_counts
        total_vol += frame_vol
    n_frames = len(per_frame_results)

    mean_vol  = total_vol / n_frames
    r_mid     = 0.5 * (bins[:-1] + bins[1:])
    shell_vol = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)
    gr        = (2 * mean_vol) / (n_proteins * (n_proteins - 1)) * total_counts / (shell_vol * n_frames)

    return {'r': r_mid, 'gr': gr, 'n_proteins': n_proteins, 'mean_vol': mean_vol}
