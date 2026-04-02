import pytest
from martinisoup.data_structures import BindingState, ResidenceRegistry


class TestBindingState:

    def test_initial_state(self):
        s = BindingState()
        assert s.bound is False
        assert s.start is None
        assert s.durations == []

    def test_moltype_id_stored(self):
        s = BindingState(moltype_id=3)
        assert s.moltype_id == 3

    def test_single_event_recorded(self):
        s = BindingState()
        s.update_bound(True, 0)
        s.update_bound(False, 10)
        assert s.durations == [10]

    def test_multiple_events_recorded(self):
        s = BindingState()
        s.update_bound(True, 0)
        s.update_bound(False, 5)
        s.update_bound(True, 10)
        s.update_bound(False, 20)
        assert s.durations == [5, 10]

    def test_repeated_bound_does_not_split_event(self):
        s = BindingState()
        s.update_bound(True, 0)
        s.update_bound(True, 5)   # still bound — should not restart event
        s.update_bound(False, 10)
        assert s.durations == [10]  # event started at 0, not 5

    def test_repeated_unbound_is_harmless(self):
        s = BindingState()
        s.update_bound(False, 0)
        s.update_bound(False, 5)
        assert s.durations == []

    def test_finalize_closes_open_event(self):
        s = BindingState()
        s.update_bound(True, 0)
        s.finalize(15)
        assert s.durations == [15]
        assert s.bound is False

    def test_finalize_while_unbound_has_no_effect(self):
        s = BindingState()
        s.update_bound(True, 0)
        s.update_bound(False, 10)
        s.finalize(20)
        assert s.durations == [10]  # no extra duration added


class TestResidenceRegistry:

    def _make_two_type_registry(self):
        """Two atoms of different molecule types, one event each."""
        s0 = BindingState(moltype_id=0)
        s0.update_bound(True, 0)
        s0.update_bound(False, 10)

        s1 = BindingState(moltype_id=1)
        s1.update_bound(True, 5)
        s1.update_bound(False, 15)

        return ResidenceRegistry(
            tracker={0: s0, 1: s1},
            molnums=[0, 1],
            moltype_table={0: 'ATP', 1: 'GTP'},
            moltype_ids=[0, 1],
        )

    def test_aggregate_by_molecule(self):
        reg = self._make_two_type_registry()
        assert reg.molecule_data[0] == [10]
        assert reg.molecule_data[1] == [10]

    def test_aggregate_by_type(self):
        reg = self._make_two_type_registry()
        assert reg.type_data[0] == [10]
        assert reg.type_data[1] == [10]

    def test_get_durations_by_type_returns_named_dict(self):
        reg = self._make_two_type_registry()
        result = reg.get_durations_by_type()
        assert set(result.keys()) == {'ATP', 'GTP'}
        assert result['ATP'] == [10]
        assert result['GTP'] == [10]

    def test_get_durations_by_type_raises_without_table(self):
        s = BindingState(moltype_id=0)
        reg = ResidenceRegistry({0: s}, [0], moltype_table=None, moltype_ids=[0])
        with pytest.raises(ValueError):
            reg.get_durations_by_type()

    def test_atoms_in_same_molecule_pool_durations(self):
        """Two atoms in molecule 0 should both contribute to molecule_data[0]."""
        s0 = BindingState(moltype_id=0)
        s0.update_bound(True, 0)
        s0.update_bound(False, 5)

        s1 = BindingState(moltype_id=0)
        s1.update_bound(True, 10)
        s1.update_bound(False, 20)

        reg = ResidenceRegistry(
            tracker={0: s0, 1: s1},
            molnums=[0, 0],   # both atoms → same molecule
            moltype_table={0: 'ATP'},
            moltype_ids=[0, 0],
        )
        assert sorted(reg.molecule_data[0]) == [5, 10]

    def test_to_dict_from_dict_roundtrip(self):
        reg = self._make_two_type_registry()
        d = reg.to_dict()
        restored = ResidenceRegistry.from_dict(d)
        assert restored.molecule_data == reg.molecule_data
        assert restored.type_data == reg.type_data
        assert restored.moltype_table == reg.moltype_table
