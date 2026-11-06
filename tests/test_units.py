import argparse
import pickle

import numpy as np
import pint
import pytest
import MDAnalysis as mda
from MDAnalysis import units as mda_units

from martinisoup.units import (
    Q_, UnitSystem, get_quantity, get_units, pint_unit_name, pretty,
    quantity_arg, unit_string, ureg,
)


def _write_xtc(tmp_path):
    """Two atoms at x = 10, 20 Å in a 50 Å box, written as .gro + .xtc (nm on disk)."""
    u = mda.Universe.empty(2, trajectory=True)
    u.add_TopologyAttr('name', ['A', 'B'])
    u.add_TopologyAttr('resname', ['X'])
    u.atoms.positions = [[10.0, 0, 0], [20.0, 0, 0]]
    u.dimensions = [50.0, 50.0, 50.0, 90, 90, 90]
    gro, xtc = str(tmp_path / 't.gro'), str(tmp_path / 't.xtc')
    u.atoms.write(gro)
    with mda.Writer(xtc, n_atoms=2) as w:
        w.write(u.atoms)
    return gro, xtc


class TestMDAnalysisUnitNames:

    @pytest.mark.parametrize("name", sorted(mda_units.lengthUnit_factor))
    def test_length_names_match_mdanalysis_factors(self, name):
        # lengthUnit_factor[name] = number of `name` units per Å
        one_angstrom = Q_(1.0, 'angstrom').to(pint_unit_name(name)).magnitude
        assert one_angstrom == pytest.approx(mda_units.lengthUnit_factor[name])

    @pytest.mark.parametrize("name", sorted(mda_units.timeUnit_factor))
    def test_time_names_match_mdanalysis_factors(self, name):
        one_ps = Q_(1.0, 'picosecond').to(pint_unit_name(name)).magnitude
        assert one_ps == pytest.approx(mda_units.timeUnit_factor[name])

    def test_A_is_angstrom_not_ampere(self):
        assert pint_unit_name('A') == 'angstrom'

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError):
            pint_unit_name('furlong')


class TestUnitSystemFromUniverse:

    def test_converted_reader_reports_base_units_not_file_units(self, tmp_path):
        """The msd bug: .xtc says nm, but MDAnalysis holds Å in memory."""
        u = mda.Universe(*_write_xtc(tmp_path))
        assert u.trajectory.units['length'] == 'nm'
        assert u.atoms.positions[1, 0] == pytest.approx(20.0)
        assert UnitSystem.from_universe(u) == UnitSystem('angstrom', 'picosecond')

    def test_unconverted_reader_reports_file_units(self, tmp_path):
        u = mda.Universe(*_write_xtc(tmp_path), convert_units=False)
        assert u.atoms.positions[1, 0] == pytest.approx(2.0)
        assert UnitSystem.from_universe(u) == UnitSystem('nanometer', 'picosecond')

    def test_memory_reader_reports_base_units(self):
        u = mda.Universe.empty(1, trajectory=True)
        assert UnitSystem.from_universe(u) == UnitSystem('angstrom', 'picosecond')

    def test_magnitudes_convert_quantities_and_pass_numbers_through(self):
        nm_system = UnitSystem('nanometer', 'picosecond')
        assert nm_system.length_magnitude(Q_(5.0, 'angstrom')) == pytest.approx(0.5)
        assert nm_system.length_magnitude(0.5) == 0.5
        assert nm_system.time_magnitude(Q_(2.0, 'ns')) == pytest.approx(2000.0)
        with pytest.raises(pint.DimensionalityError):
            nm_system.length_magnitude(Q_(1.0, 'ps'))

    def test_pickles(self):
        s = UnitSystem('nanometer', 'nanosecond')
        assert pickle.loads(pickle.dumps(s)) == s


class TestQuantityArg:

    def test_bare_number_takes_default_unit(self):
        assert quantity_arg('angstrom')('5') == Q_(5.0, 'angstrom')

    def test_value_with_units_is_kept(self):
        q = quantity_arg('angstrom')('0.5 nm')
        assert q.to('angstrom').magnitude == pytest.approx(5.0)

    def test_offset_temperature(self):
        assert quantity_arg('kelvin')('27 degC').to('kelvin').magnitude == pytest.approx(300.15)

    def test_wrong_dimension_rejected(self):
        with pytest.raises(argparse.ArgumentTypeError):
            quantity_arg('angstrom')('5 ps')

    def test_unparseable_rejected(self):
        with pytest.raises(argparse.ArgumentTypeError):
            quantity_arg('angstrom')('five')

    def test_works_inside_argparse(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--cutoff', type=quantity_arg('angstrom'), default=Q_(5.0, 'angstrom'))
        assert parser.parse_args(['--cutoff', '0.6 nm']).cutoff.to('angstrom').magnitude == pytest.approx(6.0)
        assert parser.parse_args([]).cutoff == Q_(5.0, 'angstrom')


class TestSavedData:

    def test_unit_strings_round_trip_through_pint(self):
        for unit in (ureg.angstrom ** 2, ureg.centimeter ** 2 / ureg.second, ureg.dimensionless):
            assert ureg.Unit(unit_string(unit)) == unit

    def test_pretty(self):
        assert pretty('angstrom ** 2') == 'Å²'
        assert pretty('dimensionless') == ''

    def test_get_quantity(self):
        data = {'r': [1.0, 2.0], 'units': {'r': 'nanometer'}}
        assert np.allclose(get_quantity(data, 'r').to('angstrom').magnitude, [10.0, 20.0])

    def test_legacy_data_warns_and_uses_fallback(self):
        with pytest.warns(UserWarning, match="no recorded units"):
            units = get_units({'r': [1.0]}, legacy={'r': 'angstrom'})
        assert units == {'r': 'angstrom'}

    def test_missing_units_without_fallback_raises(self):
        with pytest.raises(KeyError):
            get_units({'r': [1.0]})
