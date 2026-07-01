import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader

from martinisoup.parallel import chunk_frame_indices, map_trajectory_parallel


def _make_universe(n_frames, n_atoms=6, box_length=100.0, seed=0):
    rng = np.random.default_rng(seed)
    u = mda.Universe.empty(
        n_atoms, n_residues=n_atoms, n_segments=1,
        atom_resindex=list(range(n_atoms)), residue_segindex=[0] * n_atoms,
        trajectory=True,
    )
    u.add_TopologyAttr('name', ['CA'] * n_atoms)
    u.add_TopologyAttr('resname', ['ALA'] * n_atoms)
    u.add_TopologyAttr('resid', list(range(1, n_atoms + 1)))

    coords = rng.uniform(0, box_length, size=(n_frames, n_atoms, 3))
    dims = np.tile([box_length] * 3 + [90.0, 90.0, 90.0], (n_frames, 1))
    u.trajectory = MemoryReader(coords, order='fac', dimensions=dims, dt=1.0)
    return u


def _write_temp_files(u, tmp_path):
    gro = str(tmp_path / 'top.gro')
    xtc = str(tmp_path / 'traj.xtc')
    u.trajectory[0]
    u.atoms.write(gro)
    with mda.Writer(xtc, n_atoms=u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            w.write(u.atoms)
    return gro, xtc


def _setup_all(u, sel):
    return u.select_atoms(sel)


def _per_frame_sum(u, atoms, ts):
    return ts.frame, ts.time, float(atoms.positions.sum())


class TestChunkFrameIndices:

    def test_exact_multiple_of_chunk_size(self):
        chunks = chunk_frame_indices(0, 10, 1, chunk_size=5)
        assert chunks == [list(range(0, 5)), list(range(5, 10))]

    def test_remainder_chunk(self):
        chunks = chunk_frame_indices(0, 12, 1, chunk_size=5)
        assert [len(c) for c in chunks] == [5, 5, 2]

    def test_step_is_applied_before_chunking(self):
        chunks = chunk_frame_indices(0, 20, 4, chunk_size=2)
        flat = [i for c in chunks for i in c]
        assert flat == list(range(0, 20, 4))

    def test_empty_range_gives_no_chunks(self):
        assert chunk_frame_indices(5, 5, 1, chunk_size=10) == []

    def test_chunk_size_larger_than_range_gives_one_chunk(self):
        chunks = chunk_frame_indices(0, 3, 1, chunk_size=100)
        assert chunks == [[0, 1, 2]]


class TestMapTrajectoryParallel:

    def test_matches_serial_iteration_order_and_values(self, tmp_path):
        """
        Results must come back in ascending frame order with uneven chunk
        boundaries relative to the frame range, since several analyses
        (e.g. residence_tracker) depend on strict ordering for correctness.
        """
        u = _make_universe(53, n_atoms=6)
        gro, xtc = _write_temp_files(u, tmp_path)

        u_direct = mda.Universe(gro, xtc)
        atoms_direct = u_direct.select_atoms('all')
        serial_results = [
            _per_frame_sum(u_direct, atoms_direct, ts)
            for ts in u_direct.trajectory[3:47:2]
        ]

        parallel_results = map_trajectory_parallel(
            gro, xtc, _setup_all, _per_frame_sum,
            setup_args=('all',), start=3, stop=47, step=2,
            n_workers=3, chunk_size=7, desc="test",
        )

        assert [r[0] for r in parallel_results] == [r[0] for r in serial_results]
        assert [r[1] for r in parallel_results] == [r[1] for r in serial_results]
        assert np.allclose([r[2] for r in parallel_results], [r[2] for r in serial_results])

    def test_covers_whole_trajectory_when_stop_is_none(self, tmp_path):
        u = _make_universe(10, n_atoms=4)
        gro, xtc = _write_temp_files(u, tmp_path)

        results = map_trajectory_parallel(
            gro, xtc, _setup_all, _per_frame_sum,
            setup_args=('all',), n_workers=2, chunk_size=3, desc="test",
        )

        assert [r[0] for r in results] == list(range(10))

    def test_single_worker_single_chunk(self, tmp_path):
        u = _make_universe(4, n_atoms=4)
        gro, xtc = _write_temp_files(u, tmp_path)

        results = map_trajectory_parallel(
            gro, xtc, _setup_all, _per_frame_sum,
            setup_args=('all',), n_workers=1, chunk_size=100, desc="test",
        )

        assert [r[0] for r in results] == [0, 1, 2, 3]
