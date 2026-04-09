from __future__ import annotations

import argparse
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import PowerLawModel

from martinisoup.database import load_metabolite_classes, DATABASE_URL

COLORS = {
    'Ions':              '#337AEA',
    'Nucleotides':       '#944CE5',
    'Cofactors':         '#B7C61E',
    'Amino_Acids':       '#F9B714',
    'Other_Metabolites': '#CC7FDB',
    'Lipids':            '#05938E',
    'Carbohydrates':     '#3FB760',
}

DEFAULT_WEIGHTS_END = {
    'Ions':              5e2,
    'Nucleotides':       1e4,
    'Cofactors':         1e4,
    'Amino_Acids':       1e1,
    'Carbohydrates':     1e5,
    'Lipids':            1e3,
    'Other_Metabolites': 1e4,
}

DEFAULT_COLOR = '#888888'


def load_datasets(files: list[str], summarised: bool = True) -> list[dict]:
    datasets = []
    for f in files:
        with open(f, 'rb') as fh:
            data = pickle.load(fh)
        residences = data['residences']
        if not summarised:
            residences = {
                resname: np.unique(durations, return_counts=True)
                for resname, durations in residences.items()
            }
        datasets.append(residences)
    return datasets


def build_histograms(datasets: list[dict],
                     bins: np.ndarray,
                     metabolite_classes: dict | None = None) -> dict:
    """Aggregate per-replica residence histograms.

    If metabolite_classes is provided, histograms are grouped and averaged by
    metabolite class: {class_name: {replica_idx: {'counts': ..., 'edges': ...}}}

    Otherwise, each resname gets its own histogram:
    {resname: {replica_idx: {'counts': ..., 'edges': ...}}}
    """
    hists = defaultdict(dict)

    for idx, replica in enumerate(datasets):
        groups = defaultdict(list)

        for resname, unique_lifetimes in replica.items():
            if metabolite_classes is not None:
                if resname not in metabolite_classes:
                    continue
                key = metabolite_classes[resname]
            else:
                key = resname

            values, counts = unique_lifetimes
            hist = np.histogram(values / 1000, bins=bins, weights=counts)
            groups[key].append(hist)

        for key, vals in groups.items():
            hists[key][idx] = {
                'counts': np.stack([v[0] for v in vals]),
                'edges':  np.stack([v[1] for v in vals]),
            }

    return dict(hists)


def fit_and_plot(hists: dict, weights_end: dict, output_dir: Path,
                 grouped: bool = True) -> dict:
    """Fit a power law to each histogram group and produce plots.

    When grouped=True, keys are metabolite class names and COLORS are applied.
    When grouped=False, keys are resnames and a default colour is used.
    """
    keys = [k for k in COLORS if k in hists] if grouped else list(hists.keys())
    results = {}
    fit_data = {}  # store (x, counts, counts_std, best_fit) per key for plotting

    for key in keys:
        data = hists[key]

        # shape: (n_replicas, n_members_in_group, n_bins)
        all_counts = np.stack([r['counts'] for r in data.values()])
        edges = np.stack([r['edges'] for r in data.values()]).mean(axis=0).mean(axis=0)

        counts = all_counts.mean(axis=0).mean(axis=0)
        counts_std = all_counts.reshape(-1, all_counts.shape[-1]).std(axis=0)

        x = (edges[1:] + edges[:-1]) / 2
        x0 = x / 1000

        w_end = weights_end.get(key, 1e3)  # necessary because of log binning
        fit_weights = np.logspace(np.log10(1e0), np.log10(w_end), num=counts.shape[0])

        mod = PowerLawModel()
        params = mod.guess(counts, x=x0)
        result = mod.fit(counts, params=params, x=x0, weights=fit_weights)

        results[key] = [result.params['exponent'].value,
                        result.params['exponent'].stderr]
        fit_data[key] = (x, counts, counts_std, result.best_fit)

    # --- grid of per-group fit plots ---
    ncols = 3
    nrows = int(np.ceil(len(keys) / ncols))
    fig_fits, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                  figsize=(5 * ncols, 4 * nrows),
                                  sharey=True)
    axes_flat = np.array(axes).flatten()

    for ax, key in zip(axes_flat, keys):
        x, counts, counts_std, best_fit = fit_data[key]
        exp, err = results[key]
        color = COLORS.get(key, DEFAULT_COLOR)
        ax.errorbar(x, counts,
                    yerr=counts_std, ls='',
                    marker='o', markersize=10,
                    markeredgewidth=1, markeredgecolor='#262626',
                    c=color, label='data')
        ax.plot(x, best_fit, c='#262626',
                lw=3, ls='--', label='fit', zorder=10)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(key.replace('_', ' '))
        ax.set_xlabel("Lifetime (ns)")
        ax.set_ylabel("Count")
        ax.legend(title=f"exp = {exp:.2f} ± {err:.2f}")

    for ax in axes_flat[len(keys):]:
        ax.set_visible(False)

    fig_fits.tight_layout()
    fig_fits.savefig(output_dir / "residence_fits.png", bbox_inches='tight')
    plt.close(fig_fits)

    # --- bar chart of exponents ---
    sorted_results = dict(sorted(results.items(),
                                 key=lambda item: item[1][0],
                                 reverse=True))

    fig_bar, ax_bar = plt.subplots()
    ax_bar.bar(
        ['\n'.join(k.split('_')) for k in sorted_results],
        [v[0] for v in sorted_results.values()],
        color=[COLORS.get(k, DEFAULT_COLOR) for k in sorted_results],
        yerr=[v[1] for v in sorted_results.values()],
    )
    ax_bar.set_ylabel("Power law exponent")
    fig_bar.savefig(output_dir / "residence_exponents.png", bbox_inches='tight')
    plt.close(fig_bar)

    return results


def save_results(results: dict, output_dir: Path,
                 grouped: bool = True, command: str = '') -> None:
    csv_file = output_dir / "residence_exponents.csv"
    key_label = "class" if grouped else "resname"
    with open(csv_file, 'w') as fh:
        fh.write(f"{key_label},exponent,exponent_err\n")
        for key, (exp, err) in results.items():
            fh.write(f"{key},{exp},{err}\n")
    print(f"Results written to {csv_file}")

    pkl_file = output_dir / "residence_exponents.pkl"
    with open(pkl_file, 'wb') as fh:
        pickle.dump({"command": command, "results": results}, fh)
    print(f"Results written to {pkl_file}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit residence-time data using a power law model. "
                    "Multiple input files are treated as replicas and averaged. "
                    "If a database is provided, results are grouped by metabolite class; "
                    "otherwise each molecule type is fitted independently."
    )
    parser.add_argument(
        "files", nargs='+',
        help="Pickle file(s) from `martinisoup residence-times`. Multiple files are averaged as replicas."
    )
    parser.add_argument(
        "--database", nargs='?', const=DATABASE_URL, default=None,
        help="Group results by metabolite class using the M3-Metabolome database. "
             "Use --database alone to fetch the remote default, or supply a local CSV path."
    )
    parser.add_argument(
        "--bins-start", type=float, default=0.0,
        help="Start of log-spaced bin range in ns (default: 0, i.e. 10^0 = 1)."
    )
    parser.add_argument(
        "--bins-stop", type=float, default=np.log10(500),
        help="End of log-spaced bin range in ns as log10 value (default: log10(500))."
    )
    parser.add_argument(
        "--bins-n", type=int, default=25,
        help="Number of bins (default: 25)."
    )
    parser.add_argument(
        "--weights", default=None,
        help="Path to a JSON file mapping class names to fit weight upper bounds. "
             "Overrides built-in defaults."
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
        "--unsummarised", action="store_true",
        help="Input files are the raw (unsummarised) output of `martinisoup residence-times`. "
             "By default, summarised output (produced with --summary) is expected."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    command = ' '.join(sys.argv)

    if args.style:
        plt.style.use(args.style)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    weights_end = DEFAULT_WEIGHTS_END.copy()
    if args.weights:
        with open(args.weights) as fh:
            weights_end.update(json.load(fh))

    bins = np.logspace(args.bins_start, args.bins_stop, args.bins_n)

    metabolite_classes = None
    if args.database is not None:
        metabolite_classes = load_metabolite_classes(args.database)

    grouped = metabolite_classes is not None

    datasets = load_datasets(args.files, summarised=not args.unsummarised)
    hists = build_histograms(datasets, bins, metabolite_classes)
    results = fit_and_plot(hists, weights_end, output_dir, grouped)
    save_results(results, output_dir, grouped, command)


if __name__ == "__main__":
    main()
