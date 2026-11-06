from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import LinearModel

from martinisoup.database import load_metabolite_classes, DATABASE_URL
from martinisoup.units import Q_, get_units, pretty, quantity_arg, unit_string, ureg

# Units the fitter works in. load_datasets converts every input file to these
# from the units recorded in it, so replicas are always consistent.
FIT_LENGTH = 'angstrom'
FIT_TIME = 'picosecond'
# Units diffusion coefficients are reported in.
D_UNIT = 'centimeter ** 2 / second'

_LENGTH_SQUARED_KEYS = ('residue_timeseries', 'residue_std', 'msd')
_TIME_KEYS = ('time', 'dt', 'lagtimes')
_LENGTH_KEYS = ('mean_box_lengths',)

# `martinisoup msd` has always written MDAnalysis internal units (Å, ps), so
# files from before units were recorded are in these. (Their
# 'trajectory_units' entry describes the trajectory file, not the data, and
# must not be used.) These files also use the old key layout, in which
# 'lagtimes' held the all-particle MSD (see _upgrade_old_layout).
_LEGACY_UNITS = {
    'residue_timeseries': 'angstrom ** 2',
    'residue_std': 'angstrom ** 2',
    'lagtimes': 'angstrom ** 2',
    'time': 'picosecond',
    'dt': 'picosecond',
    'mean_box_lengths': 'angstrom',
    'dimensions': 'dimensionless',
}


def _upgrade_old_layout(data: dict, units: dict) -> None:
    """
    Convert the old MSD key layout to the current one, in place.

    Older ``martinisoup msd`` output stored the all-particle MSD under
    ``'lagtimes'`` and had no ``'msd'`` key. Current output stores it under
    ``'msd'`` and puts the actual lag times (``index × dt``) in
    ``'lagtimes'``. The old layout is recognised by the absence of ``'msd'``
    together with ``'lagtimes'`` having units of length².
    """
    if 'msd' in data or 'lagtimes' not in data:
        return
    if ureg.Unit(units['lagtimes']).dimensionality != ureg.Unit('meter ** 2').dimensionality:
        return
    data['msd'] = data.pop('lagtimes')
    units['msd'] = units.pop('lagtimes')
    data['lagtimes'] = np.arange(len(data['msd'])) * data['dt']
    units['lagtimes'] = units['dt']


def load_datasets(files: list[str]) -> list[dict]:
    """
    Load MSD pickles, converting every quantity to the fitting units.

    MSDs end up in FIT_LENGTH² (Å²), times in FIT_TIME (ps) and box lengths
    in FIT_LENGTH (Å), whatever units each file recorded; the converted
    units are stored in each dataset's ``'units'``. Files without recorded
    units are assumed to be in Å and ps, with a warning. Files in the old key
    layout (all-particle MSD under ``'lagtimes'``) are converted to the
    current one, so every dataset has ``'msd'`` (length²) and ``'lagtimes'``
    (time).
    """
    targets = {
        **{k: f'{FIT_LENGTH} ** 2' for k in _LENGTH_SQUARED_KEYS},
        **{k: FIT_TIME for k in _TIME_KEYS},
        **{k: FIT_LENGTH for k in _LENGTH_KEYS},
    }
    datasets = []
    for f in files:
        with open(f, 'rb') as fh:
            data = pickle.load(fh)
        units = dict(get_units(data, legacy=_LEGACY_UNITS, source=f))
        _upgrade_old_layout(data, units)
        for key, target in targets.items():
            if key in data:
                data[key] = Q_(np.asarray(data[key]), units[key]).to(target).magnitude
        data['units'] = {**units, **{k: unit_string(t) for k, t in targets.items()}}
        datasets.append(data)
    return datasets


def build_lagtimes(data: dict) -> np.ndarray:
    """
    Lag times of the MSD points, in FIT_TIME for datasets from load_datasets.

    Uses the saved ``'lagtimes'``; for a dataset that has not been through
    load_datasets and is in the old layout, computes ``index × dt``.
    """
    if 'msd' in data:
        return np.asarray(data['lagtimes'])
    return np.arange(len(data['lagtimes'])) * data['dt']


def yeh_hummer_correction(D, box_lengths, temperature, eta):
    """Apply the Yeh-Hummer finite-size correction to a diffusion coefficient.

    D_corr = D_PBC + xi * k_B * T / (6 * pi * eta * L)

    where L = (a * b * c)^(1/3) is the effective box length.

    Each argument may be a pint Quantity in any compatible unit. Plain
    numbers are taken to be in the units below.

    Parameters
    ----------
    D : float or pint.Quantity
        Raw (PBC) diffusion coefficient (plain number: cm²/s).
    box_lengths : array-like or pint.Quantity
        Box edge lengths [a, b, c] (plain numbers: Å).
    temperature : float or pint.Quantity
        Simulation temperature (plain number: K).
    eta : float or pint.Quantity
        Dynamic viscosity (plain number: mPa·s).

    Returns
    -------
    float or pint.Quantity
        Corrected diffusion coefficient: a Quantity if ``D`` was one,
        otherwise a plain number in cm²/s.
    """
    def as_q(value, unit):
        if isinstance(value, ureg.Quantity):
            return value
        return Q_(value, unit)

    if isinstance(box_lengths, ureg.Quantity):
        box_lengths = box_lengths.to('angstrom').magnitude

    D_q = as_q(D, D_UNIT)
    a, b, c = as_q(np.asarray(box_lengths, dtype=float), 'angstrom')
    L = (a * b * c) ** (1 / 3)
    # .to('kelvin') first, so offset units such as degC become absolute.
    T = as_q(temperature, 'kelvin').to('kelvin')
    xi = 2.837298
    correction = (xi * ureg.boltzmann_constant * T
                  / (6 * np.pi * as_q(eta, 'millipascal * second') * L))
    D_corr = D_q + correction.to(D_q.units)
    return D_corr if isinstance(D, ureg.Quantity) else D_corr.to(D_UNIT).magnitude


def average_replicas(datasets: list[dict]) -> dict:
    """Return a single dataset whose timeseries and std are averaged across replicas.

    For each residue the mean MSD is the arithmetic mean over replicas.
    The combined uncertainty is the propagated error of that mean::

        combined_std = sqrt(sum(std_i ** 2)) / n_replicas

    Datasets must share units, which is guaranteed for datasets returned by
    :func:`load_datasets`.
    """
    if len(datasets) == 1:
        return datasets[0]

    ref = datasets[0]
    n = len(datasets)

    # Build per-resname lookup for each dataset so ordering doesn't matter
    def resname_map(ds):
        return {name: (ts, st)
                for name, ts, st in zip(ds['resnames'],
                                        ds['residue_timeseries'],
                                        ds['residue_std'])}

    maps = [resname_map(ds) for ds in datasets]
    resnames = ref['resnames']

    averaged_timeseries = []
    averaged_std = []

    for name in resnames:
        ts = np.array([m[name][0] for m in maps])
        st = np.array([m[name][1] for m in maps])

        averaged_timeseries.append(ts.mean(axis=0))
        averaged_std.append(np.sqrt((st ** 2).sum(axis=0)) / n)

    return {
        'lagtimes': ref['lagtimes'],
        'msd': (np.mean([ds['msd'] for ds in datasets], axis=0)
                if all('msd' in ds for ds in datasets) else None),
        'dt': ref['dt'],
        'dimensions': ref.get('dimensions', 3),
        'resnames': resnames,
        'residue_timeseries': averaged_timeseries,
        'residue_std': averaged_std,
        'units': ref.get('units'),
    }


def fit_and_plot(data: dict, lagtimes: np.ndarray,
                 cut_start: int, cut_end: int,
                 output_dir: Path,
                 box_lengths: np.ndarray | None = None,
                 temperature=None,
                 eta=0.7103,  # mPa·s (Martini 3 water); a Quantity is also accepted
                 metabolite_classes: dict | None = None) -> dict:
    """Fit MSD curves and return results.

    If metabolite_classes is provided, results are grouped::

        {class: {resname: array([D, D_err])}}                        (no correction)
        {class: {resname: array([D_PBC, D_err, D_corr, D_err])}}    (with correction)

    Otherwise, results are flat with the same array shapes.

    The Yeh-Hummer correction is applied when both *temperature* and *box_lengths*
    are provided (see :func:`yeh_hummer_correction` for the units of plain
    numbers; pint Quantities are also accepted).

    MSDs and lag times are taken to be in FIT_LENGTH² against FIT_TIME (Å²
    against ps), which is what :func:`load_datasets` produces from any input
    units. Diffusion coefficients are returned as plain numbers in D_UNIT
    (cm²/s).
    """
    results = {}

    apply_correction = temperature is not None and box_lengths is not None

    # The fitted slope is in FIT_LENGTH² / FIT_TIME; this converts it to D_UNIT.
    slope_unit = ureg.Unit(f'{FIT_LENGTH} ** 2 / {FIT_TIME}')
    unit_factor = Q_(1.0, slope_unit).to(D_UNIT).magnitude
    d_label = pretty(D_UNIT)

    # Einstein relation: MSD = 2 * d * D * t
    einstein_factor = 2 * data.get('dimensions', 3)

    for timeseries, std, resname in zip(data['residue_timeseries'],
                                        data['residue_std'],
                                        data['resnames']):

        x_fit = lagtimes[cut_start:cut_end]
        y_fit = timeseries[cut_start:cut_end]
        y_err = std[cut_start:cut_end]

        mod = LinearModel()
        pars = mod.guess(y_fit, x=x_fit)
        result = mod.fit(y_fit, params=pars, x=x_fit, weights=1 / y_err)

        D = result.params['slope'].value / einstein_factor * unit_factor
        D_err = result.params['slope'].stderr / einstein_factor * unit_factor

        if apply_correction:
            D_corr = yeh_hummer_correction(D, box_lengths, temperature, eta)
            entry = np.array([D, D_err, D_corr, D_err])
            legend_title = (f"D_PBC = {D:.2e} ± {D_err:.2e}\n"
                            f"D_corr = {D_corr:.2e} ± {D_err:.2e} {d_label}")
        else:
            entry = np.array([D, D_err])
            legend_title = f"D = {D:.2e} ± {D_err:.2e} {d_label}"

        if metabolite_classes is not None:
            cl = metabolite_classes.get(resname, 'Unknown')
            results.setdefault(cl, {})[resname] = entry
        else:
            results[resname] = entry

        fig, ax = plt.subplots()
        ax.errorbar(lagtimes, timeseries, yerr=std,
                    ls='', marker='.', c='#93E1D1',
                    markeredgecolor='#36221A', markeredgewidth=1,
                    label='data')
        ax.plot(x_fit, result.best_fit,
                c='#3F615A', label='fit', lw=3, zorder=10)
        ax.set_title(resname, fontsize=20)
        ax.set_xlabel(f"Time ({pretty(FIT_TIME)})")
        ax.set_ylabel(f"MSD ({pretty(f'{FIT_LENGTH} ** 2')})")
        ax.legend(loc='upper left', title=legend_title)
        if metabolite_classes is not None:
            fig.savefig(output_dir / f"{cl}_{resname}.png", bbox_inches='tight')
        else:
            fig.savefig(output_dir / f"{resname}.png", bbox_inches='tight')
        plt.close(fig)

    return results


def save_results(results: dict, output_dir: Path, command: str = '') -> None:
    grouped = isinstance(next(iter(results.values())), dict)

    # Detect whether Yeh-Hummer correction was applied (4-element arrays)
    if grouped:
        _sample = next(iter(next(iter(results.values())).values()))
    else:
        _sample = next(iter(results.values()))
    corrected = len(_sample) == 4

    u = pretty(D_UNIT)
    csv_file = output_dir / "diffusion_coefficients.csv"
    with open(csv_file, 'w') as fh:
        if grouped:
            if corrected:
                fh.write(f"class,resname,D_PBC ({u}),D_err ({u}),D_corr ({u})\n")
                for cl, metabolites in results.items():
                    for resname, (D, D_err, D_corr, _) in metabolites.items():
                        fh.write(f"{cl},{resname},{D},{D_err},{D_corr}\n")
            else:
                fh.write(f"class,resname,D ({u}),D_err ({u})\n")
                for cl, metabolites in results.items():
                    for resname, (D, D_err) in metabolites.items():
                        fh.write(f"{cl},{resname},{D},{D_err}\n")
        else:
            if corrected:
                fh.write(f"resname,D_PBC ({u}),D_err ({u}),D_corr ({u})\n")
                for resname, (D, D_err, D_corr, _) in results.items():
                    fh.write(f"{resname},{D},{D_err},{D_corr}\n")
            else:
                fh.write(f"resname,D ({u}),D_err ({u})\n")
                for resname, (D, D_err) in results.items():
                    fh.write(f"{resname},{D},{D_err}\n")
    print(f"Results written to {csv_file}")

    pkl_file = output_dir / "diffusion_coefficients.pkl"
    with open(pkl_file, 'wb') as fh:
        pickle.dump({"command": command,
                     "results": results,
                     "units": {"results": unit_string(ureg.Unit(D_UNIT))}}, fh)
    print(f"Results written to {pkl_file}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit MSD data to extract diffusion coefficients. "
                    "Multiple input files are treated as replicas and averaged."
    )
    parser.add_argument(
        "files", nargs='+',
        help="Pickle file(s) containing MSD data. Multiple files are averaged as replicas."
    )
    parser.add_argument(
        "--cut-start", type=int, default=10,
        help="Start index of the fitting window (default: 10)."
    )
    parser.add_argument(
        "--cut-end", type=int, default=50,
        help="End index of the fitting window (default: 50)."
    )
    parser.add_argument(
        "--output-dir", default=".",
        help="Directory for output plots and results CSV (default: current directory)."
    )
    parser.add_argument(
        "--style", default=None,
        help="Path to a matplotlib style file."
    )
    parser.add_argument(
        "--database", nargs='?', const=DATABASE_URL, default=None,
        help="Group results by metabolite class using the M3-Metabolome database. "
             "Use --database alone to fetch the remote default, or supply a local CSV path."
    )
    parser.add_argument(
        "--temperature", type=quantity_arg('kelvin'), default=None,
        help="Simulation temperature. When provided, applies the Yeh-Hummer "
             "finite-size correction using the box dimensions stored in the MSD data. "
             "A bare number is in K; a value with units, e.g. '37 degC', is converted."
    )
    parser.add_argument(
        "--eta", type=quantity_arg('millipascal * second'), default=Q_(0.7103, 'millipascal * second'),
        help="Dynamic viscosity for the Yeh-Hummer correction "
             "(default: 0.7103 mPa·s, Martini 3 water). A bare number is in mPa·s; "
             "a value with units, e.g. '0.00071 Pa*s', is converted."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    command = ' '.join(sys.argv)

    if args.style:
        plt.style.use(args.style)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metabolite_classes = None
    if args.database is not None:
        metabolite_classes = load_metabolite_classes(args.database)

    datasets = load_datasets(args.files)

    box_lengths = None
    if args.temperature is not None:
        # load_datasets has converted every file's box lengths to FIT_LENGTH
        box_lengths = Q_(np.mean([ds['mean_box_lengths'] for ds in datasets], axis=0),
                         FIT_LENGTH)

    data = average_replicas(datasets)
    lagtimes = build_lagtimes(datasets[0])

    results = fit_and_plot(data, lagtimes, args.cut_start, args.cut_end, output_dir,
                           box_lengths=box_lengths,
                           temperature=args.temperature,
                           eta=args.eta,
                           metabolite_classes=metabolite_classes)
    save_results(results, output_dir, command)


if __name__ == "__main__":
    main()
