"""
End-to-end tests: run each ``martinisoup`` command on files on disk, as a user
would, and check the files it writes (including the recorded units).
"""

import pickle

import numpy as np
import pytest

from tests.integration.conftest import (
    D_EXPECTED_CM2_S, DT_PS, N_ATP, N_FRAMES, N_GLC, N_PROT,
)

pytestmark = pytest.mark.integration


def _load(path):
    with open(path, 'rb') as fh:
        return pickle.load(fh)


def _strip_command(data):
    return {k: v for k, v in data.items() if k != 'command'}


# --------------------------------------------------------------------------- #
# CLI dispatcher
# --------------------------------------------------------------------------- #

class TestDispatcher:

    def test_no_command_prints_usage_and_fails(self, run_cli, capsys):
        with pytest.raises(SystemExit) as exc:
            run_cli()
        assert exc.value.code == 1
        assert 'Usage: martinisoup' in capsys.readouterr().out

    def test_help_lists_commands(self, run_cli, capsys):
        with pytest.raises(SystemExit) as exc:
            run_cli('--help')
        assert exc.value.code == 0
        out = capsys.readouterr().out
        for command in ('binding-frequency', 'residence-times', 'msd', 'msd-fitter'):
            assert command in out

    def test_unknown_command_fails(self, run_cli, capsys):
        with pytest.raises(SystemExit) as exc:
            run_cli('not-a-command')
        assert exc.value.code == 1
        assert 'Unknown command' in capsys.readouterr().out

    def test_value_with_wrong_units_rejected(self, run_cli, system, tmp_path):
        with pytest.raises(SystemExit) as exc:
            run_cli('protein-rdf', system.topology, system.trajectory,
                    '--r-max', '5 ps', '--output', tmp_path / 'x.pkl')
        assert exc.value.code == 2  # argparse usage error


# --------------------------------------------------------------------------- #
# Analysis commands
# --------------------------------------------------------------------------- #

class TestProteinRdf:

    def test_serial_parallel_and_units(self, run_cli, system, tmp_path):
        common = (system.topology, system.trajectory, '--n-bins', 30)
        run_cli('protein-rdf', *common, '--r-max', '30', '--output', tmp_path / 's.pkl')
        run_cli('protein-rdf', *common, '--r-max', '3 nm', '--parallel',
                '--n_workers', 2, '--chunk_size', 7, '--output', tmp_path / 'p.pkl')
        serial, parallel = _load(tmp_path / 's.pkl'), _load(tmp_path / 'p.pkl')

        assert serial['n_proteins'] == N_PROT
        assert np.allclose(parallel['r'], serial['r'])
        assert np.allclose(parallel['gr'], serial['gr'])
        assert serial['r'].max() < 30.0
        assert serial['units'] == {
            'r': 'angstrom', 'gr': 'dimensionless',
            'n_proteins': 'dimensionless', 'mean_vol': 'angstrom ** 3',
        }


class TestProteinClustering:

    def test_serial_parallel_and_units(self, run_cli, system, tmp_path):
        common = (system.topology, system.trajectory, '--start', 5, '--stop', 50, '--step', 3)
        run_cli('protein-clustering', *common, '--r-max', '6', '--output', tmp_path / 's.pkl')
        run_cli('protein-clustering', *common, '--r-max', '0.6 nm', '--parallel',
                '--n_workers', 2, '--chunk_size', 4, '--output', tmp_path / 'p.pkl')
        serial, parallel = _load(tmp_path / 's.pkl'), _load(tmp_path / 'p.pkl')

        assert [f['frame'] for f in serial['frames']] == list(range(5, 50, 3))
        assert [f['links'] for f in parallel['frames']] == [f['links'] for f in serial['frames']]
        assert len(serial['nodes']) == N_PROT
        assert serial['r_max'] == pytest.approx(6.0) == parallel['r_max']
        assert serial['units'] == {
            'r_max': 'angstrom', 'frames.frame': 'dimensionless', 'frames.time': 'picosecond',
        }


class TestMetabolitePartitioning:

    def test_serial_parallel_and_units(self, run_cli, system, tmp_path):
        common = (system.topology, system.trajectory)
        run_cli('metabolite-partitioning', *common, '--r-max', '5', '--contact-cutoff', '6',
                '--output', tmp_path / 's.pkl')
        run_cli('metabolite-partitioning', *common, '--r-max', '0.5 nm',
                '--contact-cutoff', '600 pm', '--parallel', '--n_workers', 2,
                '--chunk_size', 9, '--output', tmp_path / 'p.pkl')
        serial, parallel = _load(tmp_path / 's.pkl'), _load(tmp_path / 'p.pkl')

        assert len(serial['results']) == N_FRAMES
        assert parallel['results'] == serial['results']
        for frame in serial['results']:
            assert set(frame['fractions']) == {'ATP', 'GLC'}
            for fractions in frame['fractions'].values():
                assert sum(fractions.values()) == pytest.approx(1.0)
        assert parallel['contact_cutoff'] == pytest.approx(6.0)
        assert serial['units']['results.time'] == 'picosecond'
        assert serial['units']['contact_cutoff'] == 'angstrom'


class TestResidenceTimes:

    def test_serial_parallel_frames_and_units(self, run_cli, system, tmp_path):
        common = (system.topology, system.trajectory, '--cutoff')
        run_cli('residence-times', *common, '6', '--output', tmp_path / 's.pkl',
                '--summary', tmp_path / 's_sum.pkl')
        run_cli('residence-times', *common, '0.6 nm', '--parallel', '--n_workers', 2,
                '--chunk_size', 7, '--output', tmp_path / 'p.pkl')
        run_cli('residence-times', *common, '6', '--frames', '--output', tmp_path / 'f.pkl')
        serial, parallel, frames = (_load(tmp_path / n) for n in ('s.pkl', 'p.pkl', 'f.pkl'))

        assert set(serial['residences']) == {'ATP', 'GLC'}
        assert sum(len(v) for v in serial['residences'].values()) > 0
        assert parallel['residences'] == serial['residences']
        assert serial['units'] == {'residences': 'picosecond', 'dt': 'picosecond',
                                   'cutoff': 'angstrom'}
        assert serial['stamp'] == 'time' and serial['dt'] == pytest.approx(DT_PS)

        # --frames: the same events, counted in frames (one frame = dt).
        assert frames['stamp'] == 'frame'
        assert frames['units']['residences'] == 'dimensionless'
        for moltype, durations in serial['residences'].items():
            assert np.allclose(np.array(frames['residences'][moltype]) * DT_PS, durations)

        summary = _load(tmp_path / 's_sum.pkl')
        values, counts = summary['residences']['ATP']
        assert sum(counts) == len(serial['residences']['ATP'])
        assert summary['units']['residences.values'] == 'picosecond'


@pytest.mark.xfail(
    strict=True, raises=TypeError,
    reason="Known bug (CLAUDE.md): binding-frequency's index_protein_segments "
           "only handles .tpr-style 'seg_<n>_<name>' segment names and crashes "
           "on other topologies.",
)
class TestBindingFrequency:

    def test_runs_and_records_units(self, run_cli, system, tmp_path):
        run_cli('binding-frequency', system.topology, system.trajectory,
                '--cutoff', '0.6 nm', '--output', tmp_path / 'b.pkl')
        out = _load(tmp_path / 'b.pkl')
        assert out['cutoff'] == pytest.approx(6.0)
        assert out['units']['cutoff'] == 'angstrom'
        assert out['n_proteins'] == {'PROT': N_PROT}
        assert out['n_metabolites'] == {'ATP': N_ATP, 'GLC': N_GLC}


# --------------------------------------------------------------------------- #
# Pipelines: analysis command -> fitter
# --------------------------------------------------------------------------- #

class TestMsdPipeline:

    def test_msd_then_fitter_recovers_known_D(self, run_cli, system, tmp_path):
        run_cli('msd', system.topology, system.trajectory,
                '--start', 0, '--stop', N_FRAMES, '--output', tmp_path / 'msd.pkl')
        msd = _load(tmp_path / 'msd.pkl')
        assert np.allclose(msd['lagtimes'], np.arange(N_FRAMES) * DT_PS)
        assert msd['units']['msd'] == 'angstrom ** 2'
        assert msd['mean_box_lengths'] == pytest.approx([60.0] * 3)  # Å, not nm

        out = tmp_path / 'fit'
        run_cli('msd-fitter', tmp_path / 'msd.pkl', '--cut-start', 5, '--cut-end', 40,
                '--temperature', '300 K', '--output-dir', out)
        fit = _load(out / 'diffusion_coefficients.pkl')

        assert fit['units'] == {'results': 'centimeter ** 2 / second'}
        for resname in ('ATP', 'GLC'):
            D, D_err, D_corr, _ = fit['results'][resname]
            # A finite random walk scatters around the exact value; 25% is
            # loose enough to be robust and tight enough to catch unit errors
            # (which are factors of 10 or more).
            assert D == pytest.approx(D_EXPECTED_CM2_S, rel=0.25)
            assert D_corr > D
            assert (out / f'{resname}.png').exists()
        header = (out / 'diffusion_coefficients.csv').read_text().splitlines()[0]
        assert header == 'resname,D_PBC (cm²/s),D_err (cm²/s),D_corr (cm²/s)'


class TestResidencePipeline:

    def test_residence_times_then_fitter(self, run_cli, system, tmp_path):
        for name, extra in (('time', ()), ('frames', ('--frames',))):
            run_cli('residence-times', system.topology, system.trajectory, '--cutoff', '8',
                    *extra, '--output', tmp_path / f'{name}.pkl',
                    '--summary', tmp_path / f'{name}_sum.pkl')

        # Lifetimes here run from 10 ps to ~1.5 ns; these bins are all populated
        # (the power-law fit cannot handle empty bins).
        bins = ('--bins-start', np.log10(0.01), '--bins-stop', np.log10(1.5), '--bins-n', 6)
        for name in ('time', 'frames'):
            run_cli('residence-fitter', tmp_path / f'{name}_sum.pkl', *bins,
                    '--output-dir', tmp_path / f'fit_{name}')
        fit_time = _load(tmp_path / 'fit_time' / 'residence_exponents.pkl')
        fit_frames = _load(tmp_path / 'fit_frames' / 'residence_exponents.pkl')

        assert fit_time['units'] == {'results': 'dimensionless'}
        assert set(fit_time['results']) == {'ATP', 'GLC'}
        # Frame-counted input is converted with dt, so the fit is identical.
        for key, (exponent, _) in fit_time['results'].items():
            assert np.isfinite(exponent)
            assert fit_frames['results'][key][0] == pytest.approx(exponent)
        assert (tmp_path / 'fit_time' / 'residence_fits.png').exists()
