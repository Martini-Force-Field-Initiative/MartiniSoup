"""
Units across the analysis classes: recorded output units, Quantity
parameters, and Universes whose data is not in Å/ps.
"""

import numpy as np
import pint
import pytest
import MDAnalysis as mda

from martinisoup.contact import ContactCounts
from martinisoup.metabolite_partitioning import MetabolitePartitioning
from martinisoup.protein_clustering import ProteinContactNetwork
from martinisoup.protein_rdf import ProteinRDF
from martinisoup.residence_tracker import ResidenceTracker
from martinisoup.units import Q_, UnitSystem

from tests.test_contact import _all_protein_types, _make_multiframe_universe
from tests.test_metabolite_partitioning import _make_random_universe, _selections
from tests.test_protein_clustering import _make_protein_universe as _make_network_universe
from tests.test_protein_rdf import _make_protein_universe as _make_rdf_universe
from tests.test_residence_tracker import (
    METAB_SEL, PROT_SEL, _BOUND_POS, _UNBOUND_POS, _make_universe, _write_temp_files,
)


class TestRecordedUnitsAndQuantityParameters:

    def test_protein_rdf(self):
        u = _make_rdf_universe(8)
        plain = ProteinRDF(u.atoms, r_max=40.0, n_bins=20).run()
        quantity = ProteinRDF(u.atoms, r_max=Q_(4.0, 'nm'), n_bins=20).run()

        assert np.allclose(quantity.results.r, plain.results.r)
        assert np.allclose(quantity.results.gr, plain.results.gr)
        assert plain.results.units == {
            'counts': 'dimensionless', 'volumes': 'angstrom ** 3',
            'r': 'angstrom', 'gr': 'dimensionless',
            'n_proteins': 'dimensionless', 'mean_vol': 'angstrom ** 3',
        }

    def test_contact_counts(self):
        u = _make_multiframe_universe(10, n_prot=3, n_met=4)
        proteins, metabolites = u.select_atoms('resname ALA'), u.select_atoms('resname ATP')
        types = _all_protein_types(proteins)
        plain = ContactCounts(proteins, metabolites, types, 20.0).run()
        quantity = ContactCounts(proteins, metabolites, types, Q_(2.0, 'nm')).run()

        assert quantity.results.counts == plain.results.counts
        assert quantity.results.cutoff_radius == pytest.approx(20.0)
        assert plain.results.units == {'counts': 'dimensionless', 'cutoff_radius': 'angstrom'}

    def test_protein_contact_network(self):
        u = _make_network_universe(5)
        plain = ProteinContactNetwork(u.atoms, r_max=6.0).run()
        quantity = ProteinContactNetwork(u.atoms, r_max=Q_(0.6, 'nm')).run()

        assert [f['links'] for f in quantity.results.frames] == \
            [f['links'] for f in plain.results.frames]
        assert plain.results.units == {
            'frames.frame': 'dimensionless', 'frames.time': 'picosecond', 'r_max': 'angstrom',
        }

    def test_metabolite_partitioning(self):
        u = _make_random_universe(5)
        plain = MetabolitePartitioning(*_selections(u), r_max=5.0, contact_cutoff=4.0).run()
        quantity = MetabolitePartitioning(
            *_selections(u), r_max=Q_(0.5, 'nm'), contact_cutoff=Q_(400, 'pm')).run()

        assert quantity.results.frames == plain.results.frames
        assert quantity.results.contact_cutoff == pytest.approx(4.0)
        assert plain.results.units == {
            'frames.frame': 'dimensionless', 'frames.time': 'picosecond',
            'frames.fractions': 'dimensionless', 'r_max': 'angstrom',
            'contact_cutoff': 'angstrom',
        }

    def test_residence_tracker_time_and_frame_stamps(self):
        u = _make_universe([_BOUND_POS, _BOUND_POS, _UNBOUND_POS])
        met, prot = u.select_atoms(METAB_SEL), u.select_atoms(PROT_SEL)
        by_time = ResidenceTracker(met, prot, cutoff=Q_(0.5, 'nm')).run()
        by_frame = ResidenceTracker(met, prot, cutoff=5.0, use_time=False).run()

        assert by_time.results.durations == by_frame.results.durations == {'ATP': [2]}
        assert by_time.results.units == {'durations': 'picosecond', 'cutoff': 'angstrom'}
        assert by_time.results.stamp == 'time'
        assert by_frame.results.units['durations'] == 'dimensionless'
        assert by_frame.results.stamp == 'frame'

    def test_wrong_dimension_rejected(self):
        u = _make_rdf_universe(2)
        with pytest.raises(pint.DimensionalityError):
            ProteinRDF(u.atoms, r_max=Q_(4.0, 'ps'))


class TestUnconvertedUniverse:
    """
    A Universe read with convert_units=False holds the file's native units
    (nm for .xtc). Analyses must work in, and report, those units, giving the
    same physical answer as a normally read Universe.
    """

    def test_residence_tracker_in_nanometres(self, tmp_path):
        frames = [_BOUND_POS, _BOUND_POS, _UNBOUND_POS, _BOUND_POS, _UNBOUND_POS]
        gro, xtc = _write_temp_files(_make_universe(frames), tmp_path)

        converted = mda.Universe(gro, xtc)
        native = mda.Universe(gro, xtc, convert_units=False)
        assert UnitSystem.from_universe(native) == UnitSystem('nanometer', 'picosecond')

        cutoff = Q_(5.0, 'angstrom')
        rt_a = ResidenceTracker(converted.select_atoms(METAB_SEL),
                                converted.select_atoms(PROT_SEL), cutoff=cutoff).run()
        rt_nm = ResidenceTracker(native.select_atoms(METAB_SEL),
                                 native.select_atoms(PROT_SEL), cutoff=cutoff).run()

        assert rt_nm.results.durations == rt_a.results.durations
        assert rt_nm.results.cutoff == pytest.approx(0.5)
        assert rt_nm.results.units['cutoff'] == 'nanometer'
        assert rt_a.results.cutoff == pytest.approx(5.0)

    def test_contact_counts_in_nanometres(self, tmp_path):
        u = _make_multiframe_universe(12, n_prot=3, n_met=5)
        gro = str(tmp_path / 'top.gro')
        xtc = str(tmp_path / 'traj.xtc')
        u.trajectory[0]
        u.atoms.write(gro)
        with mda.Writer(xtc, n_atoms=u.atoms.n_atoms) as w:
            for ts in u.trajectory:
                w.write(u.atoms)

        results = {}
        for name, kwargs in (('A', {}), ('nm', {'convert_units': False})):
            v = mda.Universe(gro, xtc, **kwargs)
            prot, met = v.select_atoms('resname ALA'), v.select_atoms('resname ATP')
            results[name] = ContactCounts(prot, met, _all_protein_types(prot),
                                          Q_(2.0, 'nm')).run().results

        assert results['nm'].counts == results['A'].counts
        assert results['nm'].units['cutoff_radius'] == 'nanometer'
        assert results['nm'].cutoff_radius == pytest.approx(2.0)
