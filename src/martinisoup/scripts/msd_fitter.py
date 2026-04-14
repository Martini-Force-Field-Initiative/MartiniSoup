from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import LinearModel
from uncertainties.unumpy import uarray, nominal_values, std_devs
from MDAnalysis.units import convert as mda_convert

from martinisoup.database import load_metabolite_classes, DATABASE_URL


def load_datasets(files: list[str]) -> list[dict]:
    datasets = []
    for f in files:
        with open(f, 'rb') as fh:
            datasets.append(pickle.load(fh))
    return datasets


def build_lagtimes(data: dict) -> np.ndarray:
    nframes = len(data['lagtimes'])
    timestep = data['dt'] 
    return np.arange(nframes) * timestep


def yeh_hummer_correction(D_cm2_s: float, box_lengths_angstrom: np.ndarray,
                          temperature: float, eta_mPas: float) -> float:
    """Apply the Yeh-Hummer finite-size correction to a diffusion coefficient.

    D_corr = D_PBC + xi * k_B * T / (6 * pi * eta * L)

    where L = (a * b * c)^(1/3) is the effective box length.

    Parameters
    ----------
    D_cm2_s : float
        Raw (PBC) diffusion coefficient in cm²/s.
    box_lengths_angstrom : array-like
        Box edge lengths [a, b, c] in Å.
    temperature : float
        Simulation temperature in K.
    eta_mPas : float
        Dynamic viscosity in mPa·s.

    Returns
    -------
    float
        Corrected diffusion coefficient in cm²/s.
    """
    xi = 2.837298
    k_B = 1.380649e-23  # J/K
    a, b, c = box_lengths_angstrom
    L_m = (a * b * c) ** (1 / 3) * 1e-10  # Å → m
    eta_Pa_s = eta_mPas * 1e-3             # mPa·s → Pa·s
    correction_m2_s = xi * k_B * temperature / (6 * np.pi * eta_Pa_s * L_m)
    return D_cm2_s + correction_m2_s * 1e4  # m²/s → cm²/s


def average_replicas(datasets: list[dict]) -> dict:
    """Return a single dataset whose timeseries and std are averaged across replicas.

    For each residue the mean MSD is the arithmetic mean over replicas.
    The combined uncertainty is the propagated error of that mean::

        combined_std = sqrt(sum(std_i ** 2)) / n_replicas
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
        'dt': ref['dt'],
        'resnames': resnames,
        'residue_timeseries': averaged_timeseries,
        'residue_std': averaged_std,
    }


def fit_and_plot(data: dict, lagtimes: np.ndarray,
                 cut_start: int, cut_end: int,
                 output_dir: Path,
                 time_unit: str, length_unit: str,
                 box_lengths: np.ndarray | None = None,
                 temperature: float | None = None,
                 eta: float = 0.7103,
                 metabolite_classes: dict | None = None) -> dict:
    """Fit MSD curves and return results.

    If metabolite_classes is provided, results are grouped:
        {class: {resname: array([D, D_err])}}                        (no correction)
        {class: {resname: array([D_PBC, D_err, D_corr, D_err])}}    (with correction)
    Otherwise, results are flat with the same array shapes.

    The Yeh-Hummer correction is applied when both *temperature* and *box_lengths*
    are provided.
    """
    results = {}

    apply_correction = temperature is not None and box_lengths is not None

    # slope is in length_unit²/time_unit; convert to cm²/s using stored trajectory units.
    # MDAnalysis does not recognise 'cm' or 's', so normalise to Angstrom and ps first,
    # then apply exact physical constants (1 Å = 1e-8 cm, 1 ps = 1e-12 s).
    length_in_angstrom = mda_convert(1.0, length_unit, 'Angstrom')
    time_in_ps = mda_convert(1.0, time_unit, 'ps')
    unit_factor = (length_in_angstrom * 1e-8) ** 2 / (time_in_ps * 1e-12)

    for timeseries, std, resname in zip(data['residue_timeseries'],
                                        data['residue_std'],
                                        data['resnames']):

        x_fit = lagtimes[cut_start:cut_end]
        y_fit = timeseries[cut_start:cut_end]
        y_err = std[cut_start:cut_end]

        mod = LinearModel()
        pars = mod.guess(y_fit, x=x_fit)
        result = mod.fit(y_fit, params=pars, x=x_fit, weights=1 / y_err ** 2)

        D = result.params['slope'].value / 6 * unit_factor
        D_err = result.params['slope'].stderr / 6 * unit_factor

        if apply_correction:
            D_corr = yeh_hummer_correction(D, box_lengths, temperature, eta)
            entry = np.array([D, D_err, D_corr, D_err])
            legend_title = (f"D_PBC = {D:.2e} ± {D_err:.2e}\n"
                            f"D_corr = {D_corr:.2e} ± {D_err:.2e} cm²/s")
        else:
            entry = np.array([D, D_err])
            legend_title = f"D = {D:.2e} ± {D_err:.2e} cm²/s"

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
        ax.set_xlabel(f"Time ({time_unit})")
        ax.set_ylabel(f"MSD ({length_unit}²)")
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

    csv_file = output_dir / "diffusion_coefficients.csv"
    with open(csv_file, 'w') as fh:
        if grouped:
            if corrected:
                fh.write("class,resname,D_PBC (cm²/s),D_err (cm²/s),D_corr (cm²/s)\n")
                for cl, metabolites in results.items():
                    for resname, (D, D_err, D_corr, _) in metabolites.items():
                        fh.write(f"{cl},{resname},{D},{D_err},{D_corr}\n")
            else:
                fh.write("class,resname,D (cm²/s),D_err (cm²/s)\n")
                for cl, metabolites in results.items():
                    for resname, (D, D_err) in metabolites.items():
                        fh.write(f"{cl},{resname},{D},{D_err}\n")
        else:
            if corrected:
                fh.write("resname,D_PBC (cm²/s),D_err (cm²/s),D_corr (cm²/s)\n")
                for resname, (D, D_err, D_corr, _) in results.items():
                    fh.write(f"{resname},{D},{D_err},{D_corr}\n")
            else:
                fh.write("resname,D (cm²/s),D_err (cm²/s)\n")
                for resname, (D, D_err) in results.items():
                    fh.write(f"{resname},{D},{D_err}\n")
    print(f"Results written to {csv_file}")

    pkl_file = output_dir / "diffusion_coefficients.pkl"
    with open(pkl_file, 'wb') as fh:
        pickle.dump({"command": command, 
                     "results": results,
                     "units": "cm²/s"}, fh)
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
        "--temperature", type=float, default=300,
        help="Simulation temperature in K. When provided, applies the Yeh-Hummer "
             "finite-size correction using the box dimensions stored in the MSD data."
    )
    parser.add_argument(
        "--eta", type=float, default=0.7103,
        help="Dynamic viscosity in mPa·s for the Yeh-Hummer correction "
             "(default: 0.7103 mPa·s, Martini 3 water)."
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
    traj_units = datasets[0]['trajectory_units']
    time_unit = traj_units['time']
    length_unit = traj_units['length']

    box_lengths = None
    if args.temperature is not None:
        box_lengths = np.mean([ds['mean_box_lengths'] for ds in datasets], axis=0)

    data = average_replicas(datasets)
    lagtimes = build_lagtimes(datasets[0])

    results = fit_and_plot(data, lagtimes, args.cut_start, args.cut_end, output_dir,
                           time_unit, length_unit,
                           box_lengths=box_lengths,
                           temperature=args.temperature,
                           eta=args.eta,
                           metabolite_classes=metabolite_classes)
    save_results(results, output_dir, command)


if __name__ == "__main__":
    main()
