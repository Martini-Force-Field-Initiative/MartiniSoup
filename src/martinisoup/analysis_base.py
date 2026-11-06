"""
Shared base class for martinisoup analyses built on MDAnalysis's AnalysisBase.

Subclasses get the 'serial', 'multiprocessing' and 'dask' backends of
``AnalysisBase.run``. The built-in 'multiprocessing' backend is swapped for
one that dispatches exactly one computation part per pool task.

Why: MDAnalysis's ``BackendMultiprocessing`` (checked against 2.10 and 2.11) calls
``Pool.map(partial(self._compute), parts)`` with the default chunksize,
``ceil(n_parts / (4 * n_workers))``. When that exceeds 1, several parts are
sent in one batch and run on the *same* unpickled analysis object; each
``_compute`` overwrites ``self.results`` and returns ``self``, so every part
in the batch comes back holding the last part's results. Results are silently
wrong whenever ``n_parts > 4 * n_workers``. Forcing ``chunksize=1`` gives each
part its own copy of the object.
"""

import math

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.backends import BackendMultiprocessing

from .units import UnitSystem


def n_parts_for_chunk_size(trajectory, start, stop, step, chunk_size):
    """
    Number of parts so each holds roughly ``chunk_size`` frames.

    Lets the ``*_parallel`` wrappers keep their ``chunk_size`` argument while
    passing ``n_parts`` to ``AnalysisBase.run``.
    """
    n_frames = len(trajectory[start:stop:step])
    return max(1, math.ceil(n_frames / chunk_size))


class _BackendMultiprocessingOnePartPerTask(BackendMultiprocessing):
    """``BackendMultiprocessing`` with ``chunksize=1`` (see module docstring)."""

    def apply(self, func, computations):
        from multiprocessing import Pool

        with Pool(processes=self.n_workers) as pool:
            return pool.map(func, computations, chunksize=1)


class ParallelAnalysisBase(AnalysisBase):
    """
    AnalysisBase for analyses whose per-frame work is independent.

    Subclasses implement the usual ``_prepare`` / ``_single_frame`` /
    ``_conclude`` plus ``_get_aggregator``, storing everything that must
    survive a parallel run in ``self.results``.

    ``self.unit_system`` holds the length/time units of the trajectory's data
    in memory (see :class:`martinisoup.units.UnitSystem`). Subclasses convert
    unitful parameters to it with ``self.unit_system.length_magnitude(...)``,
    and record the units of their outputs in ``self.results.units`` (set in
    ``_conclude``, so it is not part of the parallel aggregation).
    """

    _analysis_algorithm_is_parallelizable = True

    def __init__(self, trajectory, **kwargs):
        super().__init__(trajectory, **kwargs)
        self.unit_system = UnitSystem.from_trajectory(trajectory)

    @classmethod
    def get_supported_backends(cls):
        return ('serial', 'multiprocessing', 'dask')

    def _configure_backend(self, backend, n_workers, unsupported_backend=False):
        if isinstance(backend, str) and backend == 'multiprocessing':
            backend = _BackendMultiprocessingOnePartPerTask(n_workers=n_workers)
            unsupported_backend = True
        return super()._configure_backend(
            backend=backend,
            n_workers=n_workers,
            unsupported_backend=unsupported_backend,
        )
