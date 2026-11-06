import numpy as np
from MDAnalysis.analysis.base import ResultsGroup
from MDAnalysis.lib.distances import self_distance_array

from .analysis_base import ParallelAnalysisBase, n_parts_for_chunk_size
from .units import unit_string


class ProteinRDF(ParallelAnalysisBase):
    """
    Protein centre-of-mass radial distribution function.

    Each protein molecule's centre of mass is computed per frame, then all
    pairwise (minimum-image) distances are histogrammed. The RDF is
    normalised by the mean box volume and the number of protein pairs.

    Parallel execution is provided by MDAnalysis's AnalysisBase (via
    :class:`martinisoup.analysis_base.ParallelAnalysisBase`): the frame
    range is split into contiguous parts, each part accumulates its own
    histogram and per-frame volumes, and the parts are merged (counts summed,
    volumes concatenated) before `_conclude` normalises once in the main
    process.

    Parameters
    ----------
    atomgroup : MDAnalysis.AtomGroup
        Protein atoms. Molecules are inferred via atomgroup.split('molecule').
    r_max : float or pint.Quantity
        Upper distance limit (default: 100.0). A plain number is in the
        Universe's length unit (Å for a Universe read with default settings).
    n_bins : int
        Number of histogram bins (default: 200).
    **kwargs
        Passed to AnalysisBase (e.g. verbose).

    Results
    -------
    All lengths are in the Universe's length unit L (Å by default); the
    unit of every result is recorded in ``results.units``.

    results.counts : np.ndarray
        Raw pair counts per bin, summed over frames.
    results.volumes : np.ndarray
        Box volume of each analysed frame (L^3).
    results.r : np.ndarray
        Bin midpoints (L).
    results.gr : np.ndarray
        g(r) values (dimensionless).
    results.n_proteins : int
        Number of protein molecules.
    results.mean_vol : float
        Mean box volume over analysed frames (L^3).
    results.units : dict
        {result name: pint unit string}.

    Examples
    --------
    >>> rdf = ProteinRDF(u.select_atoms(prot_sel), r_max=100.0, n_bins=200)
    >>> rdf.run(backend='multiprocessing', n_workers=4)
    >>> rdf.results.r, rdf.results.gr
    """

    def __init__(self, atomgroup, r_max=100.0, n_bins=200, **kwargs):
        super().__init__(atomgroup.universe.trajectory, **kwargs)
        self._molecules = atomgroup.split('molecule')
        self.n_proteins = len(self._molecules)
        r_max = self.unit_system.length_magnitude(r_max)
        self.bins = np.linspace(0, r_max, n_bins + 1)

    def _prepare(self):
        self.results.counts = np.zeros(len(self.bins) - 1, dtype=np.float64)
        self.results.volumes = np.zeros(self.n_frames, dtype=np.float64)

    def _single_frame(self):
        ts = self._ts
        coms = np.array([mol.center_of_mass(wrap=True) for mol in self._molecules])
        dists = self_distance_array(coms, box=ts.dimensions)
        self.results.counts += np.histogram(dists, bins=self.bins)[0]
        # _frame_index is local to this part in parallel runs, which is what
        # we want: each part fills its own volumes array, then they are
        # concatenated in part order.
        self.results.volumes[self._frame_index] = np.prod(ts.dimensions[:3])

    def _get_aggregator(self):
        return ResultsGroup(lookup={
            'counts': ResultsGroup.ndarray_sum,
            'volumes': ResultsGroup.ndarray_hstack,
        })

    def _conclude(self):
        bins = self.bins
        n_frames = len(self.results.volumes)
        n_proteins = self.n_proteins

        mean_vol = self.results.volumes.mean()
        shell_vol = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)

        self.results.r = 0.5 * (bins[:-1] + bins[1:])
        self.results.gr = (
            (2 * mean_vol) / (n_proteins * (n_proteins - 1))
            * self.results.counts / (shell_vol * n_frames)
        )
        self.results.n_proteins = n_proteins
        self.results.mean_vol = mean_vol

        length = self.unit_system.length_unit
        self.results.units = {
            'counts': 'dimensionless',
            'volumes': unit_string(length ** 3),
            'r': unit_string(length),
            'gr': 'dimensionless',
            'n_proteins': 'dimensionless',
            'mean_vol': unit_string(length ** 3),
        }


def _as_dict(rdf):
    return {
        'r': rdf.results.r,
        'gr': rdf.results.gr,
        'n_proteins': rdf.results.n_proteins,
        'mean_vol': rdf.results.mean_vol,
        'units': {k: rdf.results.units[k] for k in ('r', 'gr', 'n_proteins', 'mean_vol')},
    }


def compute_rdf(u, atoms, r_max=100.0, n_bins=200, start=0, stop=None, step=1):
    """
    Compute the protein centre-of-mass RDF over a trajectory slice (serial).

    Thin wrapper around :class:`ProteinRDF` run with the serial backend.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Unused beyond consistency with other analyses; the trajectory is taken
        from ``atoms.universe``.
    atoms : MDAnalysis.AtomGroup
        Protein atoms. Molecules are inferred via atoms.split('molecule').
    r_max : float or pint.Quantity
        Upper distance limit (default: 100.0). A plain number is in the
        Universe's length unit (Å for a Universe read with default settings).
    n_bins : int
        Number of histogram bins (default: 200).
    start, stop, step : int
        Frame range, as in ``u.trajectory[start:stop:step]``.

    Returns
    -------
    dict with keys:
        'r'          : np.ndarray, bin midpoints (length)
        'gr'         : np.ndarray, g(r) values (dimensionless)
        'n_proteins' : int, number of protein molecules
        'mean_vol'   : float, mean box volume over analysed frames (length^3)
        'units'      : dict, pint unit string for each of the keys above
    """
    rdf = ProteinRDF(atoms, r_max=r_max, n_bins=n_bins, verbose=True)
    rdf.run(start=start, stop=stop, step=step)
    return _as_dict(rdf)


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

    Thin wrapper around :class:`ProteinRDF` run with MDAnalysis's
    ``multiprocessing`` backend. The frame range is split into contiguous
    parts of roughly ``chunk_size`` frames, distributed over ``n_workers``
    processes.

    Parameters
    ----------
    topology, trajectory : str
        Paths to topology and trajectory files.
    prot_sel : str
        MDAnalysis selection string for protein atoms.
    r_max : float or pint.Quantity
        Upper distance limit (default: 100.0). A plain number is in the
        Universe's length unit (Å, as the Universe is read with default
        settings).
    n_bins : int
        Number of histogram bins (default: 200).
    start, stop, step : int
        Frame range, as in ``u.trajectory[start:stop:step]``.
    n_workers : int
        Number of worker processes (default: 4).
    chunk_size : int
        Approximate frames per part (default: 100).

    Returns
    -------
    dict
        Same keys as :func:`compute_rdf`.
    """
    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    n_parts = n_parts_for_chunk_size(u.trajectory, start, stop, step, chunk_size)

    rdf = ProteinRDF(u.select_atoms(prot_sel), r_max=r_max, n_bins=n_bins)
    rdf.run(
        start=start, stop=stop, step=step,
        backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
    )
    return _as_dict(rdf)
