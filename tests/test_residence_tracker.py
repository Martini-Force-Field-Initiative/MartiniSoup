import pytest
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from martinisoup.residence_tracker import (
    ResidenceTracker, track_serial, track_parallel, _moltypes, _molnums,
)

# Protein sits at (5, 5, 5) throughout.  With cutoff=5 Å:
_BOUND_POS = (7.0, 5.0, 5.0)    # 2 Å from protein → bound
_UNBOUND_POS = (25.0, 5.0, 5.0) # 20 Å from protein → unbound

METAB_SEL = 'resname ATP'
PROT_SEL  = 'resname ALA'


def _make_universe(metabolite_positions_per_frame, box_length=50.0, dt=1.0):
    """
    2-atom universe: atom 0 is a protein atom fixed at (5, 5, 5), atom 1 is
    a metabolite atom whose position varies per frame.
    """
    n_frames = len(metabolite_positions_per_frame)
    protein_pos = np.array([5.0, 5.0, 5.0])

    u = mda.Universe.empty(
        2,
        n_residues=2,
        n_segments=2,
        atom_resindex=[0, 1],
        residue_segindex=[0, 1],
        trajectory=True,
    )
    u.add_TopologyAttr('name', ['CA', 'C1'])
    u.add_TopologyAttr('resname', ['ALA', 'ATP'])
    u.add_TopologyAttr('resid', [1, 2])
    u.add_TopologyAttr('moltypes', ['ALA', 'ATP'])
    u.add_TopologyAttr('molnums', [0, 1])

    coords = np.array(
        [[protein_pos, np.array(pos)] for pos in metabolite_positions_per_frame],
        dtype=float,
    )  # shape (n_frames, 2, 3)

    dimensions = np.tile(
        [box_length, box_length, box_length, 90.0, 90.0, 90.0], (n_frames, 1)
    )
    u.trajectory = MemoryReader(coords, order='fac', dimensions=dimensions, dt=dt)
    return u


def _write_temp_files(u, tmp_path):
    """Write universe topology and trajectory to temp GRO + XTC files."""
    gro = str(tmp_path / 'top.gro')
    xtc = str(tmp_path / 'traj.xtc')
    u.trajectory[0]
    u.atoms.write(gro)
    with mda.Writer(xtc, n_atoms=u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            w.write(u.atoms)
    return gro, xtc


def _make_universe_without_mol_attrs():
    """Minimal 2-residue universe with no moltypes or molnums (GRO-like)."""
    u = mda.Universe.empty(2, n_residues=2, n_segments=2,
                           atom_resindex=[0, 1], residue_segindex=[0, 1],
                           trajectory=True)
    u.add_TopologyAttr('name', ['CA', 'C1'])
    u.add_TopologyAttr('resname', ['ALA', 'ATP'])
    u.add_TopologyAttr('resid', [1, 2])
    return u


class TestMoltypeHelpers:

    def test_moltypes_returns_moltypes_when_present(self):
        u = _make_universe([_BOUND_POS])
        ag = u.select_atoms('resname ATP')
        assert list(_moltypes(ag)) == ['ATP']

    def test_moltypes_falls_back_to_resnames(self):
        u = _make_universe_without_mol_attrs()
        ag = u.select_atoms('resname ATP')
        assert list(_moltypes(ag)) == ['ATP']

    def test_molnums_returns_molnums_when_present(self):
        u = _make_universe([_BOUND_POS])
        ag = u.select_atoms('resname ATP')
        assert list(_molnums(ag)) == [1]

    def test_molnums_falls_back_to_resindices(self):
        u = _make_universe_without_mol_attrs()
        ag = u.select_atoms('resname ATP')
        # resindex of the ATP atom is 1 (second residue)
        assert list(_molnums(ag)) == list(ag.resindices)


class TestTrackSerial:

    def test_bound_metabolite_produces_duration(self, tmp_path):
        """Metabolite always within cutoff — total duration should be non-zero."""
        u = _make_universe([_BOUND_POS, _BOUND_POS, _BOUND_POS])
        gro, xtc = _write_temp_files(u, tmp_path)
        result = track_serial(gro, xtc, METAB_SEL, PROT_SEL, cutoff=5.0, use_time=False)
        assert sum(result.get('ATP', [])) > 0

    def test_always_unbound_metabolite_has_no_durations(self, tmp_path):
        """Metabolite always outside cutoff — no binding events should be recorded."""
        u = _make_universe([_UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS])
        gro, xtc = _write_temp_files(u, tmp_path)
        result = track_serial(gro, xtc, METAB_SEL, PROT_SEL, cutoff=5.0, use_time=False)
        assert result.get('ATP', []) == []

    def test_result_keys_match_molecule_types(self, tmp_path):
        u = _make_universe([_BOUND_POS, _UNBOUND_POS])
        gro, xtc = _write_temp_files(u, tmp_path)
        result = track_serial(gro, xtc, METAB_SEL, PROT_SEL, cutoff=5.0, use_time=False)
        assert set(result.keys()) == {'ATP'}

    def test_binding_event_duration_is_correct(self, tmp_path):
        """
        Frames 0,1 bound; frames 2,3,4 unbound.
        With use_time=False, frame stamps are 0,1,2,3,4.
        Binding starts at stamp 0, ends at stamp 2 → duration 2.
        """
        frames = [_BOUND_POS, _BOUND_POS, _UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS]
        u = _make_universe(frames)
        gro, xtc = _write_temp_files(u, tmp_path)
        result = track_serial(gro, xtc, METAB_SEL, PROT_SEL, cutoff=5.0, use_time=False)
        assert result.get('ATP') == [2]

    def test_start_stop_restricts_frame_range(self, tmp_path):
        """Analyzing only unbound frames should record no binding events."""
        frames = [_BOUND_POS, _BOUND_POS, _UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS]
        u = _make_universe(frames)
        gro, xtc = _write_temp_files(u, tmp_path)
        result = track_serial(gro, xtc, METAB_SEL, PROT_SEL,
                               cutoff=5.0, use_time=False, start=2, stop=5)
        assert result.get('ATP', []) == []

    def test_step_affects_recorded_duration(self, tmp_path):
        """
        Frames: bound, unbound, unbound, unbound (stamps 0,1,2,3).
        step=1: binding ends at stamp 1 → duration 1.
        step=3: frames 0 and 3 only; binding starts at stamp 0, ends at stamp 3 → duration 3.
        """
        frames = [_BOUND_POS, _UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS]
        u = _make_universe(frames)
        gro, xtc = _write_temp_files(u, tmp_path)
        result_step1 = track_serial(gro, xtc, METAB_SEL, PROT_SEL,
                                    cutoff=5.0, use_time=False, step=1)
        result_step3 = track_serial(gro, xtc, METAB_SEL, PROT_SEL,
                                    cutoff=5.0, use_time=False, step=3)
        assert result_step1.get('ATP') == [1]
        assert result_step3.get('ATP') == [3]


class TestTrackParallel:
    """
    track_parallel splits the trajectory across worker processes for the
    geometry pass, then reassembles bound-atom sets in frame order and folds
    them through the same state machine track_serial uses. These tests
    exercise that reassembly with deliberately uneven chunk boundaries.
    """

    def test_matches_serial_with_mixed_bound_unbound_pattern(self, tmp_path):
        frames = [_BOUND_POS, _BOUND_POS, _UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS] * 10
        u = _make_universe(frames)
        gro, xtc = _write_temp_files(u, tmp_path)

        serial = track_serial(gro, xtc, METAB_SEL, PROT_SEL, cutoff=5.0, use_time=False)
        parallel = track_parallel(gro, xtc, METAB_SEL, PROT_SEL, cutoff=5.0, use_time=False,
                                  n_workers=2, chunk_size=7)  # chunk boundary falls mid-event
        assert serial == parallel

    def test_matches_serial_when_binding_spans_chunk_boundary(self, tmp_path):
        """A single long bound run straddling several chunk boundaries must
        still be recorded as one contiguous duration, not split at chunk
        edges."""
        frames = ([_UNBOUND_POS] * 3 + [_BOUND_POS] * 20 + [_UNBOUND_POS] * 3)
        u = _make_universe(frames)
        gro, xtc = _write_temp_files(u, tmp_path)

        serial = track_serial(gro, xtc, METAB_SEL, PROT_SEL, cutoff=5.0, use_time=False)
        parallel = track_parallel(gro, xtc, METAB_SEL, PROT_SEL, cutoff=5.0, use_time=False,
                                  n_workers=3, chunk_size=4)
        assert serial == parallel == {'ATP': [20]}

    def test_respects_start_stop_step(self, tmp_path):
        frames = [_BOUND_POS, _BOUND_POS, _UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS]
        u = _make_universe(frames)
        gro, xtc = _write_temp_files(u, tmp_path)

        serial = track_serial(gro, xtc, METAB_SEL, PROT_SEL,
                              cutoff=5.0, use_time=False, start=2, stop=5)
        parallel = track_parallel(gro, xtc, METAB_SEL, PROT_SEL,
                                  cutoff=5.0, use_time=False, start=2, stop=5,
                                  n_workers=2, chunk_size=2)
        assert serial == parallel == {'ATP': []}


def _make_multi_metabolite_universe(n_frames, n_met=6, box_length=30.0, seed=0):
    """
    One protein atom at the box centre plus `n_met` metabolite atoms that
    random-walk, giving many bind/unbind events of varying length, some of
    them across the periodic boundary. Two moltypes, ATP and GLC.
    """
    rng = np.random.default_rng(seed)
    n_atoms = 1 + n_met
    u = mda.Universe.empty(
        n_atoms, n_residues=n_atoms, n_segments=2,
        atom_resindex=list(range(n_atoms)),
        residue_segindex=[0] + [1] * n_met,
        trajectory=True,
    )
    u.add_TopologyAttr('name', ['CA'] + ['C1'] * n_met)
    met_names = ['ATP', 'GLC'] * (n_met // 2)
    u.add_TopologyAttr('resname', ['ALA'] + met_names)
    u.add_TopologyAttr('resid', list(range(1, n_atoms + 1)))
    u.add_TopologyAttr('moltypes', ['ALA'] + met_names)
    u.add_TopologyAttr('molnums', list(range(n_atoms)))

    centre = np.full(3, box_length / 2)
    start = centre + rng.uniform(-6, 6, size=(n_met, 3))
    steps = rng.normal(0, 1.0, size=(n_frames, n_met, 3))
    met = start + np.cumsum(steps, axis=0)  # unwrapped; may leave the box
    prot = np.broadcast_to(centre, (n_frames, 1, 3))
    coords = np.concatenate([prot, met], axis=1)
    dims = np.tile([box_length] * 3 + [90.0, 90.0, 90.0], (n_frames, 1))
    u.trajectory = MemoryReader(coords, order='fac', dimensions=dims, dt=2.0)
    return u


class TestResidenceTracker:

    def test_known_event_durations_with_time_and_frame_stamps(self):
        frames = [_BOUND_POS, _BOUND_POS, _UNBOUND_POS, _BOUND_POS, _UNBOUND_POS]
        u = _make_universe(frames, dt=2.5)
        met, prot = u.select_atoms(METAB_SEL), u.select_atoms(PROT_SEL)

        by_frame = ResidenceTracker(met, prot, cutoff=5.0, use_time=False).run()
        assert by_frame.results.durations == {'ATP': [2, 1]}
        assert 'frame_contacts' not in by_frame.results

        by_time = ResidenceTracker(met, prot, cutoff=5.0, use_time=True).run()
        assert by_time.results.durations == {'ATP': [5.0, 2.5]}

    def test_does_not_add_transformations_to_universe(self):
        u = _make_universe([_BOUND_POS, _UNBOUND_POS])
        ResidenceTracker(u.select_atoms(METAB_SEL), u.select_atoms(PROT_SEL),
                         cutoff=5.0).run()
        assert len(u.trajectory.transformations) == 0

    def test_contact_across_periodic_boundary(self):
        # Protein at (5,5,5); metabolite at (48,5,5) is 7 Å away through the
        # boundary of the 50 Å box, and at (-1,5,5) it is outside the box.
        u = _make_universe([(48.0, 5.0, 5.0), (-1.0, 5.0, 5.0), _UNBOUND_POS])
        met, prot = u.select_atoms(METAB_SEL), u.select_atoms(PROT_SEL)
        rt = ResidenceTracker(met, prot, cutoff=8.0, use_time=False).run()
        assert rt.results.durations == {'ATP': [2]}

    @pytest.mark.parametrize("n_workers,n_parts", [(2, 2), (3, 7), (2, 24)])
    def test_multiprocessing_matches_serial(self, n_workers, n_parts):
        """
        Durations, including events spanning part boundaries, must match
        serial. n_parts > 4 * n_workers is the regression case for the stock
        MDAnalysis multiprocessing backend (see martinisoup.analysis_base).
        """
        u = _make_multi_metabolite_universe(101, seed=3)
        met = u.select_atoms('resname ATP GLC')
        prot = u.select_atoms(PROT_SEL)

        serial = ResidenceTracker(met, prot, cutoff=5.0).run(start=3, stop=99, step=2)
        parallel = ResidenceTracker(met, prot, cutoff=5.0).run(
            start=3, stop=99, step=2,
            backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
        )

        # Guard against a vacuous comparison: several events of each type,
        # and some longer than 2 sampled frames (step * dt * 2 = 8.0), which
        # therefore cross part boundaries in the 24-part (2 frames/part) run.
        assert all(len(serial.results.durations[k]) >= 3 for k in ('ATP', 'GLC'))
        assert max(max(v) for v in serial.results.durations.values()) > 8.0

        assert parallel.results.durations == serial.results.durations
