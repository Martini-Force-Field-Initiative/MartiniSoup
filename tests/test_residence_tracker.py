import pytest
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from martinisoup.residence_tracker import BindingEventTracker

# Protein sits at (5, 5, 5) throughout.  With cutoff=5 Å:
_BOUND_POS = (7.0, 5.0, 5.0)    # 2 Å from protein → bound
_UNBOUND_POS = (25.0, 5.0, 5.0) # 20 Å from protein → unbound


def _make_universe(metabolite_positions_per_frame, box_length=50.0):
    """
    2-atom universe: atom 0 is a protein atom fixed at (5, 5, 5), atom 1 is
    a metabolite atom whose position varies per frame.

    Parameters
    ----------
    metabolite_positions_per_frame : list of (x, y, z)
    box_length : float
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
    u.add_TopologyAttr('resid', [1, 1])
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


class TestBindingEventTracker:

    def test_bound_metabolite_produces_duration(self):
        """Metabolite always within cutoff — total duration should be non-zero."""
        u = _make_universe([_BOUND_POS, _BOUND_POS, _BOUND_POS])
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = BindingEventTracker(u, metabolites, proteins,
                                     cutoff=5.0, use_time=False).track()
        assert sum(result.get('ATP', [])) > 0

    def test_always_unbound_metabolite_has_no_durations(self):
        """Metabolite always outside cutoff — no binding events should be recorded."""
        u = _make_universe([_UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS])
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = BindingEventTracker(u, metabolites, proteins,
                                     cutoff=5.0, use_time=False).track()
        assert result.get('ATP', []) == []

    def test_result_keys_match_molecule_types(self):
        u = _make_universe([_BOUND_POS, _UNBOUND_POS])
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = BindingEventTracker(u, metabolites, proteins,
                                     cutoff=5.0, use_time=False).track()
        assert set(result.keys()) == {'ATP'}

    def test_binding_event_duration_is_correct(self):
        """
        Frames 0,1 bound; frames 2,3,4 unbound.
        With use_time=False, frame stamps are 0,1,2,3,4.
        Binding starts at stamp 0, ends at stamp 2 → duration 2.
        """
        frames = [_BOUND_POS, _BOUND_POS, _UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS]
        u = _make_universe(frames)
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = BindingEventTracker(u, metabolites, proteins,
                                     cutoff=5.0, use_time=False).track()
        assert result.get('ATP') == [2]

    def test_start_stop_restricts_frame_range(self):
        """Analyzing only unbound frames should record no binding events."""
        frames = [_BOUND_POS, _BOUND_POS, _UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS]
        u = _make_universe(frames)
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = BindingEventTracker(u, metabolites, proteins,
                                     cutoff=5.0, use_time=False,
                                     start=2, stop=5).track()
        assert result.get('ATP', []) == []

    def test_step_affects_recorded_duration(self):
        """
        Frames: bound, unbound, unbound, unbound (stamps 0,1,2,3).
        step=1: binding ends at stamp 1 → duration 1.
        step=3: frames 0 and 3 only; binding starts at stamp 0, ends at stamp 3 → duration 3.
        """
        frames = [_BOUND_POS, _UNBOUND_POS, _UNBOUND_POS, _UNBOUND_POS]
        u = _make_universe(frames)
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')

        result_step1 = BindingEventTracker(u, metabolites, proteins,
                                           cutoff=5.0, use_time=False, step=1).track()
        result_step3 = BindingEventTracker(u, metabolites, proteins,
                                           cutoff=5.0, use_time=False, step=3).track()

        assert result_step1.get('ATP') == [1]
        assert result_step3.get('ATP') == [3]
