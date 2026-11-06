import numpy as np
import pytest
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader

from martinisoup.metabolite_partitioning import (
    MetabolitePartitioning,
    analyse_trajectory,
    analyse_trajectory_parallel,
)


def _make_universe(coords, met_resnames, n_prot, box_length):
    """
    Universe of `n_prot` single-atom protein residues followed by single-atom
    metabolite residues, each residue its own molecule.

    coords : (n_frames, n_atoms, 3) array
    """
    n_frames, n_atoms, _ = coords.shape
    n_met = n_atoms - n_prot
    u = mda.Universe.empty(
        n_atoms, n_residues=n_atoms, n_segments=2,
        atom_resindex=np.arange(n_atoms),
        residue_segindex=[0] * n_prot + [1] * n_met,
        trajectory=True,
    )
    u.add_TopologyAttr('name', ['BB'] * n_prot + ['C1'] * n_met)
    u.add_TopologyAttr('resname', ['ALA'] * n_prot + list(met_resnames))
    u.add_TopologyAttr('resid', list(range(1, n_atoms + 1)))
    u.add_TopologyAttr('molnums', np.arange(n_atoms))  # per-residue attr
    dims = np.tile([box_length] * 3 + [90.0, 90.0, 90.0], (n_frames, 1))
    u.trajectory = MemoryReader(coords, order='fac', dimensions=dims, dt=2.0)
    return u


def _make_random_universe(n_frames, n_prot=4, n_met=24, box_length=30.0, seed=0):
    """Random positions in a box small enough to give a mix of all three states."""
    rng = np.random.default_rng(seed)
    coords = rng.uniform(0, box_length, size=(n_frames, n_prot + n_met, 3))
    resnames = ['ATP', 'GLC'] * (n_met // 2)
    return _make_universe(coords, resnames, n_prot, box_length)


def _selections(u):
    return u.select_atoms('resname ATP GLC'), u.select_atoms('resname ALA')


class TestMetabolitePartitioning:

    def test_known_states_single_frame(self):
        """
        Hand-placed metabolites with known states:

        - ATP at 2 Å from the protein               -> protein_adsorbed
        - ATP far from everything                   -> soluble
        - two GLC 3 Å apart                         -> clustered
        - two GLC 3 Å apart across the x boundary   -> clustered
        """
        positions = np.array([[
            [50.0, 50.0, 50.0],   # protein
            [52.0, 50.0, 50.0],   # ATP, adsorbed
            [80.0, 80.0, 80.0],   # ATP, soluble
            [10.0, 10.0, 10.0],   # GLC, clustered
            [13.0, 10.0, 10.0],   # GLC, clustered
            [1.0, 30.0, 30.0],    # GLC, clustered via PBC
            [98.0, 30.0, 30.0],   # GLC, clustered via PBC
        ]])
        u = _make_universe(positions, ['ATP', 'ATP', 'GLC', 'GLC', 'GLC', 'GLC'],
                           n_prot=1, box_length=100.0)
        metabolites, proteins = _selections(u)

        mp = MetabolitePartitioning(metabolites, proteins, r_max=5.0,
                                    contact_cutoff=5.0).run()

        (frame,) = mp.results.frames
        assert frame['frame'] == 0
        assert frame['fractions'] == {
            'ATP': {'protein_adsorbed': 0.5, 'soluble': 0.5, 'clustered': 0.0},
            'GLC': {'protein_adsorbed': 0.0, 'soluble': 0.0, 'clustered': 1.0},
        }

    def test_frames_and_times_follow_slice(self):
        u = _make_random_universe(10)
        mp = MetabolitePartitioning(*_selections(u)).run(start=1, stop=9, step=3)

        assert [f['frame'] for f in mp.results.frames] == [1, 4, 7]
        assert [f['time'] for f in mp.results.frames] == [2.0, 8.0, 14.0]

    @pytest.mark.parametrize("n_workers,n_parts", [(2, 2), (3, 7), (2, 24)])
    def test_multiprocessing_matches_serial(self, n_workers, n_parts):
        """
        Frames must come back complete and in order for uneven part
        boundaries. n_parts > 4 * n_workers is the regression case for the
        stock MDAnalysis multiprocessing backend (see martinisoup.analysis_base).
        """
        u = _make_random_universe(53)
        metabolites, proteins = _selections(u)
        serial = MetabolitePartitioning(metabolites, proteins).run(start=3, stop=50, step=2)
        parallel = MetabolitePartitioning(metabolites, proteins).run(
            start=3, stop=50, step=2,
            backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
        )

        # Guard against a vacuous comparison: all three states occur.
        states = {k for f in serial.results.frames
                  for fr in f['fractions'].values() for k, v in fr.items() if v > 0}
        assert states == {'protein_adsorbed', 'clustered', 'soluble'}

        assert [f['frame'] for f in parallel.results.frames] == list(range(3, 50, 2))
        assert parallel.results.frames == serial.results.frames


class TestWrappers:

    def test_analyse_trajectory_returns_frames(self):
        u = _make_random_universe(6)
        metabolites, proteins = _selections(u)
        frames = analyse_trajectory(u, metabolites, proteins, start=1, stop=5)

        assert [f['frame'] for f in frames] == [1, 2, 3, 4]
        assert set(frames[0]) == {'frame', 'time', 'fractions'}
        assert set(frames[0]['fractions']) == {'ATP', 'GLC'}

    def test_parallel_wrapper_matches_serial(self, monkeypatch):
        # The wrapper takes file paths, but molnums needs a .tpr, so hand it
        # the in-memory universe instead.
        u = _make_random_universe(40)
        monkeypatch.setattr(mda, 'Universe', lambda *args, **kwargs: u)
        metabolites, proteins = _selections(u)

        serial = analyse_trajectory(u, metabolites, proteins, start=2, stop=38)
        parallel = analyse_trajectory_parallel(
            'top', 'traj', 'resname ATP GLC', 'resname ALA', start=2, stop=38,
            n_workers=2, chunk_size=4,
        )
        assert parallel == serial
