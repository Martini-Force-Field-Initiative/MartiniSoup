"""Unit handling in the downstream fitters (msd-fitter, residence-fitter)."""

import pickle
import sys

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from martinisoup.scripts import msd_fitter, residence_fitter
from martinisoup.units import Q_


def _dump(path, obj):
    with open(path, 'wb') as fh:
        pickle.dump(obj, fh)
    return str(path)


# --------------------------------------------------------------------------- #
# msd-fitter
# --------------------------------------------------------------------------- #

D_TRUE = 1.0e-6        # cm²/s
D_A2_PS = 1.0e-2       # the same in Å²/ps (1 Å²/ps = 1e-4 cm²/s)
N_LAG, DT_PS = 60, 10.0


def _msd_data(length='angstrom', time='picosecond'):
    """
    3D Einstein MSD for D_TRUE (plus small fixed noise), expressed in the
    given units, in the format written by `martinisoup msd`.
    """
    t_ps = np.arange(N_LAG) * DT_PS
    # Small fixed noise: a perfect line leaves lmfit unable to estimate stderr.
    noise = np.random.default_rng(0).normal(0, 0.05, N_LAG)
    msd_a2 = 6 * D_A2_PS * t_ps + noise
    l2 = Q_(1.0, 'angstrom ** 2').to(f'{length} ** 2').magnitude
    tf = Q_(1.0, 'picosecond').to(time).magnitude
    lf = Q_(1.0, 'angstrom').to(length).magnitude
    return {
        'resnames': ['ATP', 'GLC'],
        'residue_timeseries': np.vstack([msd_a2, msd_a2]) * l2,
        'residue_std': np.full((2, N_LAG), 0.5) * l2,
        'time': t_ps * tf,
        'lagtimes': t_ps * tf,
        'msd': msd_a2 * l2,
        'dimensions': 3,
        'dt': DT_PS * tf,
        'mean_box_lengths': np.array([100.0, 100.0, 100.0]) * lf,
        'units': {
            'residue_timeseries': f'{length} ** 2', 'residue_std': f'{length} ** 2',
            'time': time, 'lagtimes': time, 'msd': f'{length} ** 2',
            'dimensions': 'dimensionless',
            'dt': time, 'mean_box_lengths': length,
        },
    }


def _old_layout(data):
    """The key layout older `martinisoup msd` wrote: all-particle MSD under 'lagtimes'."""
    data = dict(data, units=dict(data['units']))
    data['lagtimes'] = data.pop('msd')
    data['units']['lagtimes'] = data['units'].pop('msd')
    return data


def _old_yeh_hummer(D_cm2_s, box_a, T, eta_mPas):
    """The hand-coded formula msd_fitter used before pint, for comparison."""
    a, b, c = box_a
    L_m = (a * b * c) ** (1 / 3) * 1e-10
    corr = 2.837298 * 1.380649e-23 * T / (6 * np.pi * eta_mPas * 1e-3 * L_m)
    return D_cm2_s + corr * 1e4


class TestYehHummer:

    def test_plain_numbers_match_previous_formula(self):
        box = np.array([150.0, 150.0, 160.0])
        expected = _old_yeh_hummer(2e-6, box, 310.0, 0.7103)
        assert msd_fitter.yeh_hummer_correction(2e-6, box, 310.0, 0.7103) == \
            pytest.approx(expected, rel=1e-9)

    def test_quantities_in_other_units(self):
        box = np.array([150.0, 150.0, 160.0])
        expected = _old_yeh_hummer(2e-6, box, 310.0, 0.7103)
        result = msd_fitter.yeh_hummer_correction(
            Q_(2e-10, 'm ** 2 / s'), Q_(box / 10, 'nm'), Q_(36.85, 'degC'), Q_(7.103e-4, 'Pa * s'),
        )
        assert result.to('cm ** 2 / s').magnitude == pytest.approx(expected, rel=1e-9)


class TestMsdLoadDatasets:

    def test_converts_to_fitting_units(self, tmp_path):
        f = _dump(tmp_path / 'nm.pkl', _msd_data('nanometer', 'nanosecond'))
        (ds,) = msd_fitter.load_datasets([f])
        ref = _msd_data()
        for key in ('residue_timeseries', 'residue_std', 'msd', 'lagtimes', 'time', 'dt',
                    'mean_box_lengths'):
            assert np.allclose(ds[key], ref[key]), key
        assert ds['units']['residue_timeseries'] == 'angstrom ** 2'
        assert ds['units']['dt'] == 'picosecond'

    def test_legacy_file_assumes_angstrom_ps_and_ignores_trajectory_units(self, tmp_path):
        legacy = _old_layout(_msd_data())
        del legacy['units']
        legacy['trajectory_units'] = {'length': 'nm', 'time': 'ps'}  # describes the file, not the data
        f = _dump(tmp_path / 'legacy.pkl', legacy)
        with pytest.warns(UserWarning, match="no recorded units"):
            (ds,) = msd_fitter.load_datasets([f])
        assert np.allclose(ds['residue_timeseries'], _msd_data()['residue_timeseries'])
        # Old layout upgraded: MSD moved to 'msd', real lag times computed from dt.
        assert np.allclose(ds['msd'], _msd_data()['msd'])
        assert np.allclose(ds['lagtimes'], np.arange(N_LAG) * DT_PS)

    def test_old_layout_with_recorded_units_is_upgraded(self, tmp_path):
        f = _dump(tmp_path / 'old.pkl', _old_layout(_msd_data('nanometer', 'nanosecond')))
        (ds,) = msd_fitter.load_datasets([f])
        assert np.allclose(ds['msd'], _msd_data()['msd'])
        assert np.allclose(ds['lagtimes'], np.arange(N_LAG) * DT_PS)
        assert ds['units']['lagtimes'] == 'picosecond'
        assert ds['units']['msd'] == 'angstrom ** 2'


class TestMsdFit:

    def test_recovers_D_from_replicas_in_different_units(self, tmp_path):
        files = [
            _dump(tmp_path / 'a.pkl', _msd_data('angstrom', 'picosecond')),
            _dump(tmp_path / 'b.pkl', _msd_data('nanometer', 'nanosecond')),
        ]
        datasets = msd_fitter.load_datasets(files)
        data = msd_fitter.average_replicas(datasets)
        lagtimes = msd_fitter.build_lagtimes(datasets[0])
        results = msd_fitter.fit_and_plot(data, lagtimes, 10, 50, tmp_path)

        # Both replicas hold the same physical data, so the average equals
        # either one; D must come out the same as from the Å/ps file alone.
        (single,) = msd_fitter.load_datasets(files[:1])
        D_single = msd_fitter.fit_and_plot(
            single, msd_fitter.build_lagtimes(single), 10, 50, tmp_path)['ATP'][0]
        for resname in ('ATP', 'GLC'):
            assert results[resname][0] == pytest.approx(D_single, rel=1e-9)
            assert results[resname][0] == pytest.approx(D_TRUE, rel=1e-2)

    def test_main_writes_units(self, tmp_path, monkeypatch):
        f = _dump(tmp_path / 'msd.pkl', _msd_data('nanometer', 'picosecond'))
        out = tmp_path / 'out'
        monkeypatch.setattr(sys, 'argv', [
            'msd-fitter', f, '--output-dir', str(out), '--temperature', '27 degC',
        ])
        msd_fitter.main()

        with open(out / 'diffusion_coefficients.pkl', 'rb') as fh:
            saved = pickle.load(fh)
        assert saved['units'] == {'results': 'centimeter ** 2 / second'}
        D, D_err, D_corr, _ = saved['results']['ATP']
        assert D == pytest.approx(D_TRUE, rel=1e-2)
        # Box lengths were saved in nm (10 nm) and must be read as 100 Å.
        assert D_corr == pytest.approx(
            _old_yeh_hummer(D, [100.0, 100.0, 100.0], 300.15, 0.7103), rel=1e-9)

        header = (out / 'diffusion_coefficients.csv').read_text().splitlines()[0]
        assert header == 'resname,D_PBC (cm²/s),D_err (cm²/s),D_corr (cm²/s)'


# --------------------------------------------------------------------------- #
# residence-fitter
# --------------------------------------------------------------------------- #

def _residence_summary(unit, dt_ps=20.0):
    """Summary pickle; with unit='dimensionless', durations are frame counts (dt apart)."""
    values_ps = np.array([2000.0, 5000.0, 40000.0])
    if unit == 'dimensionless':
        values = values_ps / dt_ps
    else:
        values = values_ps * Q_(1.0, 'picosecond').to(unit).magnitude
    return {
        'residences': {'ATP': (values, np.array([5, 3, 1]))},
        'stamp': 'frame' if unit == 'dimensionless' else 'time',
        'dt': dt_ps,
        'cutoff': 5.0,
        'units': {'residences.values': unit, 'residences.counts': 'dimensionless',
                  'dt': 'picosecond', 'cutoff': 'angstrom'},
    }


class TestResidenceLoadDatasets:

    def test_summarised_files_converted_to_ns(self, tmp_path):
        ps = _dump(tmp_path / 'ps.pkl', _residence_summary('picosecond'))
        ns = _dump(tmp_path / 'ns.pkl', _residence_summary('nanosecond'))
        a, b = residence_fitter.load_datasets([ps, ns])
        for ds in (a, b):
            values, counts = ds['ATP']
            assert np.allclose(values, [2.0, 5.0, 40.0])
            assert list(counts) == [5, 3, 1]

    def test_unsummarised_file(self, tmp_path):
        raw = {'residences': {'ATP': [2000.0, 2000.0, 5000.0]},
               'units': {'residences': 'picosecond', 'cutoff': 'angstrom'}}
        (ds,) = residence_fitter.load_datasets([_dump(tmp_path / 'raw.pkl', raw)],
                                               summarised=False)
        values, counts = ds['ATP']
        assert np.allclose(values, [2.0, 5.0]) and list(counts) == [2, 1]

    def test_legacy_file_assumes_ps(self, tmp_path):
        legacy = _residence_summary('picosecond')
        del legacy['units']
        with pytest.warns(UserWarning, match="no recorded units"):
            (ds,) = residence_fitter.load_datasets([_dump(tmp_path / 'old.pkl', legacy)])
        assert np.allclose(ds['ATP'][0], [2.0, 5.0, 40.0])

    def test_frame_counted_durations_converted_with_dt(self, tmp_path):
        f = _dump(tmp_path / 'frames.pkl', _residence_summary('dimensionless', dt_ps=20.0))
        (ds,) = residence_fitter.load_datasets([f])
        assert np.allclose(ds['ATP'][0], [2.0, 5.0, 40.0])

    def test_frame_counted_durations_without_dt_rejected(self, tmp_path):
        data = _residence_summary('dimensionless')
        del data['dt'], data['units']['dt']
        with pytest.raises(ValueError, match="not a time"):
            residence_fitter.load_datasets([_dump(tmp_path / 'frames.pkl', data)])
