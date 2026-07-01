"""
Generic parallel map over MDAnalysis trajectory frames.

Every parallel analysis in this package has the same shape: split a
trajectory into frame-index chunks, dispatch each chunk to a worker process
that rebuilds the Universe (Universes aren't picklable) and steps through
its assigned frames, then collect per-frame results back in the main
process. This module owns that shape so individual analyses only need to
supply the per-frame science: a `setup_fn(u, *setup_args) -> context` run
once per worker, and a `per_frame_fn(u, context, ts) -> result` run once per
frame.

`setup_fn` and `per_frame_fn` must be module-level functions (not
closures/lambdas), since they are pickled to worker processes as part of the
task arguments.
"""

from concurrent.futures import ProcessPoolExecutor

from tqdm import tqdm


def chunk_frame_indices(start, stop, step, chunk_size):
    """Split a frame range into a list of index chunks of at most chunk_size frames."""
    frame_indices = list(range(start, stop, step))
    return [frame_indices[i:i + chunk_size] for i in range(0, len(frame_indices), chunk_size)]


def _run_chunk(args):
    """Worker: rebuild the Universe, run setup_fn once, then per_frame_fn per frame."""
    topology, trajectory, frame_indices, setup_fn, per_frame_fn, setup_args = args

    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    context = setup_fn(u, *setup_args)

    results = []
    for frame_idx in frame_indices:
        u.trajectory[frame_idx]
        results.append(per_frame_fn(u, context, u.trajectory.ts))
    return results


def map_trajectory_parallel(
    topology,
    trajectory,
    setup_fn,
    per_frame_fn,
    setup_args=(),
    start=0,
    stop=None,
    step=1,
    n_workers=4,
    chunk_size=100,
    desc="Processing (parallel)",
):
    """
    Compute one result per trajectory frame, in parallel, in frame order.

    Parameters
    ----------
    topology, trajectory : str
        Paths to topology and trajectory files.
    setup_fn : callable
        `setup_fn(u, *setup_args) -> context`, run once per worker process
        after rebuilding the Universe (e.g. atom selections, static metadata).
    per_frame_fn : callable
        `per_frame_fn(u, context, ts) -> result`, run once per frame.
    setup_args : tuple
        Extra arguments forwarded to `setup_fn`.
    start, stop, step : int
        Frame range, as in `u.trajectory[start:stop:step]`. `stop=None`
        means end of trajectory.
    n_workers : int
        Number of worker processes.
    chunk_size : int
        Frames per chunk dispatched to a worker.
    desc : str
        tqdm progress bar description.

    Returns
    -------
    list
        Flat list of per-frame results, in ascending frame order.
    """
    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    stop = stop if stop is not None else len(u.trajectory)

    chunks = chunk_frame_indices(start, stop, step, chunk_size)
    worker_args = [
        (topology, trajectory, chunk, setup_fn, per_frame_fn, setup_args)
        for chunk in chunks
    ]

    all_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for chunk_results in tqdm(executor.map(_run_chunk, worker_args),
                                   total=len(chunks), desc=desc):
            all_results.extend(chunk_results)
    return all_results
