import pytest
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from martinisoup.residence_tracker import track_serial, _moltypes, _molnums

# Protein sits at (5, 5, 5) throughout.  With cutoff=5 Å:
_BOUND_POS = (7.0, 5.0, 5.0)    # 2 Å from protein → bound
_UNBOUND_POS = (25.0, 5.0, 5.0) # 20 Å from protein → unbound

METAB_SEL = 'resname ATP'
PROT_SEL  = 'resname ALA'


def _make_universe(metabolite_positions_per_frame, box_length=50.0):
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
    u.trajectory = MemoryReader(coords, order='fac', dimensions=dimensions, dt=1.0)
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
