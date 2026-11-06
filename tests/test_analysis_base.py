import multiprocessing

from martinisoup import analysis_base
from martinisoup.analysis_base import _BackendMultiprocessingOnePartPerTask


def test_worker_processes_are_spawned_not_forked(monkeypatch):
    """
    Workers must be spawned whatever the platform default is. Forked workers
    (the Linux default before Python 3.14) inherit open trajectory readers and
    threads, which segfaulted in CI; see the analysis_base module docstring.
    """
    requested = []
    real_get_context = multiprocessing.get_context

    def recording_get_context(method=None):
        requested.append(method)
        return real_get_context(method)

    monkeypatch.setattr(analysis_base.multiprocessing, 'get_context', recording_get_context)

    backend = _BackendMultiprocessingOnePartPerTask(n_workers=2)
    assert backend.apply(abs, [-1, -2, -3]) == [1, 2, 3]
    assert requested == ['spawn']
