import argparse
import json
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import PowerLawModel

from martinisoup.database import load_metabolite_classes

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


def build_class_histograms(datasets: list[dict],
                           metabolite_classes: dict,
                           bins: np.ndarray) -> dict:
    """Aggregate per-replica residence histograms by metabolite class.

    Returns a nested dict: {class_name: {replica_idx: {'counts': ..., 'edges': ...}}}
    """
    hists = defaultdict(dict)

    for idx, replica in enumerate(datasets):
        survivals = defaultdict(list)

        for resname, unique_lifetimes in replica.items():
            if resname not in metabolite_classes:
                continue
            metabolite_class = metabolite_classes[resname]
            values, counts = unique_lifetimes
            hist = np.histogram(values / 1000, bins=bins, weights=counts)
            survivals[metabolite_class].append(hist)

        for cl, vals in survivals.items():
            hists[cl][idx] = {
                'counts': np.stack([v[0] for v in vals]),
                'edges':  np.stack([v[1] for v in vals]),
            }

    return dict(hists)


def fit_and_plot(hists: dict, weights_end: dict, output_dir: Path) -> dict:
    """Fit a power law to each class's averaged histogram and produce plots."""
    classes = [cl for cl in COLORS if cl in hists]
    results = {}
    fit_data = {}  # store (x, counts, counts_std, best_fit) per class for plotting

    for cl in classes:
        data = hists[cl]

        # shape: (n_replicas, n_metabolites_in_class, n_bins)
        all_counts = np.stack([r['counts'] for r in data.values()])
        edges = np.stack([r['edges'] for r in data.values()]).mean(axis=0).mean(axis=0)

        counts = all_counts.mean(axis=0).mean(axis=0)
        counts_std = all_counts.reshape(-1, all_counts.shape[-1]).std(axis=0)

        x = (edges[1:] + edges[:-1]) / 2
        x0 = x / 1000

        w_end = weights_end.get(cl, 1e3)  # necessary because of log binning
        fit_weights = np.logspace(np.log10(1e0), np.log10(w_end), num=counts.shape[0])

        mod = PowerLawModel()
        params = mod.guess(counts, x=x0)
        result = mod.fit(counts, params=params, x=x0, weights=fit_weights)

        results[cl] = [result.params['exponent'].value,
                       result.params['exponent'].stderr]
        fit_data[cl] = (x, counts, counts_std, result.best_fit)

    # --- grid of per-class fit plots ---
    ncols = 3
    nrows = int(np.ceil(len(classes) / ncols))
    fig_fits, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                  figsize=(5 * ncols, 4 * nrows),
                                  sharey=True)
    axes_flat = np.array(axes).flatten()

    for ax, cl in zip(axes_flat, classes):
        x, counts, counts_std, best_fit = fit_data[cl]
        exp, err = results[cl]
        ax.errorbar(x, counts, 
                    yerr=counts_std, ls='',
                    marker='o', markersize=10,
                    markeredgewidth=1, markeredgecolor='#262626',
                    c=COLORS[cl], label='data'
                    )
        ax.plot(x, best_fit, c='#262626', #c=COLORS[cl], 
                lw=3, ls='--', label='fit', zorder=10)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(cl.replace('_', ' '))
        ax.set_xlabel("Lifetime (ns)")
        ax.set_ylabel("Count")
        ax.legend(title=f"exp = {exp:.2f} ± {err:.2f}")

    for ax in axes_flat[len(classes):]:
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
        ['\n'.join(cl.split('_')) for cl in sorted_results],
        [v[0] for v in sorted_results.values()],
        color=[COLORS[cl] for cl in sorted_results],
        yerr=[v[1] for v in sorted_results.values()],
    )
    ax_bar.set_ylabel("Power law exponent")
    fig_bar.savefig(output_dir / "residence_exponents.png", bbox_inches='tight')
    plt.close(fig_bar)

    return results


def save_results(results: dict, output_dir: Path) -> None:
    out_file = output_dir / "residence_exponents.csv"
    with open(out_file, 'w') as fh:
        fh.write("class,exponent,exponent_err\n")
        for cl, (exp, err) in results.items():
            fh.write(f"{cl},{exp},{err}\n")
    print(f"Results written to {out_file}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit residence-time data by metabolite class using a power law model. "
                    "Multiple input files are treated as replicas and averaged."
    )
    parser.add_argument(
        "files", nargs='+',
        help="Pickle file(s) from `martinisoup residence-times`. Multiple files are averaged as replicas."
    )
    parser.add_argument(
        "--database-url", default=None,
        help="URL to metabolite class database CSV (default: M3-Metabolome repository)."
    )
    parser.add_argument(
        "--database", default=None,
        help="Path to a local metabolite class database CSV."
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

    if args.style:
        plt.style.use(args.style)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    weights_end = DEFAULT_WEIGHTS_END.copy()
    if args.weights:
        with open(args.weights) as fh:
            weights_end.update(json.load(fh))

    bins = np.logspace(args.bins_start, args.bins_stop, args.bins_n)

    datasets = load_datasets(args.files, summarised=not args.unsummarised)
    metabolite_classes = load_metabolite_classes(args.database_url, args.database)
    hists = build_class_histograms(datasets, metabolite_classes, bins)
    results = fit_and_plot(hists, weights_end, output_dir)
    save_results(results, output_dir)


if __name__ == "__main__":
    main()
