import itertools

import numpy as np
import pytest
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.lib.distances import distance_array

from martinisoup.protein_clustering import (
    ProteinContactNetwork,
    analyse_trajectory,
    analyse_trajectory_parallel,
    build_nodes,
)


def _make_protein_universe(n_frames, n_mols=6, atoms_per_mol=3, box_length=40.0, seed=0):
    """In-memory universe of small 'proteins' with molnums and moltypes."""
    rng = np.random.default_rng(seed)
    n_atoms = n_mols * atoms_per_mol
    u = mda.Universe.empty(
        n_atoms, n_residues=n_mols, n_segments=1,
        atom_resindex=np.repeat(np.arange(n_mols), atoms_per_mol),
        residue_segindex=[0] * n_mols,
        trajectory=True,
    )
    u.add_TopologyAttr('name', ['BB'] * n_atoms)
    u.add_TopologyAttr('resname', ['ALA'] * n_mols)
    u.add_TopologyAttr('resid', list(range(1, n_mols + 1)))
    u.add_TopologyAttr('molnums', np.arange(n_mols))  # per-residue attr
    u.add_TopologyAttr('moltypes', ['PA', 'PB'] * (n_mols // 2))

    # Small box relative to molecule count so each frame has a mix of
    # contacts and non-contacts, including across periodic boundaries.
    centres = rng.uniform(0, box_length, size=(n_frames, n_mols, 1, 3))
    jitter = rng.normal(0, 1.5, size=(n_frames, n_mols, atoms_per_mol, 3))
    coords = (centres + jitter).reshape(n_frames, n_atoms, 3)
    dims = np.tile([box_length] * 3 + [90.0, 90.0, 90.0], (n_frames, 1))
    u.trajectory = MemoryReader(coords, order='fac', dimensions=dims, dt=2.0)
    return u


def _reference_links(u, atoms, r_max, ts):
    """Brute-force molecule pairs with any atom pair within r_max (PBC)."""
    molnums = atoms.molnums
    links = set()
    for a, b in itertools.combinations(np.unique(molnums), 2):
        d = distance_array(atoms.positions[molnums == a],
                           atoms.positions[molnums == b], box=ts.dimensions)
        if (d <= r_max).any():
            links.add((int(a), int(b)))
    return links


def _links(frame):
    return {(l['source'], l['target']) for l in frame['links']}


class TestProteinContactNetwork:

    def test_serial_matches_brute_force(self):
        u = _make_protein_universe(10)
        atoms = u.atoms
        net = ProteinContactNetwork(atoms, r_max=6.0).run(start=1, stop=9, step=2)

        frames = net.results.frames
        assert [f['frame'] for f in frames] == [1, 3, 5, 7]
        assert [f['time'] for f in frames] == [2.0, 6.0, 10.0, 14.0]
        n_links = 0
        for f in frames:
            u.trajectory[f['frame']]
            expected = _reference_links(u, atoms, 6.0, u.trajectory.ts)
            assert _links(f) == expected
            n_links += len(expected)
        # Guard against a vacuous test: some, but not all, pairs in contact.
        assert 0 < n_links < 4 * 15

    def test_nodes_are_built_and_shared(self):
        u = _make_protein_universe(3)
        net = ProteinContactNetwork(u.atoms, r_max=6.0).run()

        assert net.results.nodes == build_nodes(u.atoms)
        assert [n['type'] for n in net.results.nodes] == ['PA', 'PB'] * 3
        assert all(f['nodes'] == net.results.nodes for f in net.results.frames)

    @pytest.mark.parametrize("n_workers,n_parts", [(2, 2), (3, 7), (2, 24)])
    def test_multiprocessing_matches_serial(self, n_workers, n_parts):
        """
        Frames must come back complete and in order for uneven part
        boundaries. n_parts > 4 * n_workers is the regression case for the
        stock MDAnalysis multiprocessing backend (see martinisoup.analysis_base).
        """
        u = _make_protein_universe(53)
        serial = ProteinContactNetwork(u.atoms, r_max=6.0).run(start=3, stop=50, step=2)
        parallel = ProteinContactNetwork(u.atoms, r_max=6.0).run(
            start=3, stop=50, step=2,
            backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
        )

        assert [f['frame'] for f in parallel.results.frames] == list(range(3, 50, 2))
        assert [f['time'] for f in parallel.results.frames] == \
            [f['time'] for f in serial.results.frames]
        assert [_links(f) for f in parallel.results.frames] == \
            [_links(f) for f in serial.results.frames]


class TestWrappers:

    def test_analyse_trajectory_returns_frames(self):
        u = _make_protein_universe(6)
        nodes = build_nodes(u.atoms)
        frames = analyse_trajectory(u, u.atoms, nodes, r_max=6.0, start=1, stop=5)

        assert [f['frame'] for f in frames] == [1, 2, 3, 4]
        assert set(frames[0]) == {
            'directed', 'multigraph', 'graph', 'nodes', 'links', 'frame', 'time',
        }

    def test_parallel_wrapper_matches_serial(self, monkeypatch):
        # The wrapper takes file paths, but molnums/moltypes need a .tpr, so
        # hand it the in-memory universe instead.
        u = _make_protein_universe(40)
        monkeypatch.setattr(mda, 'Universe', lambda *args, **kwargs: u)

        serial = analyse_trajectory(u, u.atoms, build_nodes(u.atoms), r_max=6.0,
                                    start=2, stop=38, step=1)
        parallel = analyse_trajectory_parallel(
            'top', 'traj', 'all', r_max=6.0, start=2, stop=38,
            n_workers=2, chunk_size=4,
        )
        assert [f['frame'] for f in parallel] == [f['frame'] for f in serial]
        assert [_links(f) for f in parallel] == [_links(f) for f in serial]
