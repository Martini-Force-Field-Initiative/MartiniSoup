"""
Utility functions for histograms, binning, and bootstrapping.
"""

import numpy as np


def binned_histogram(durations, nbins=20, log_bins=True):
    """
    Compute a histogram (linear or log-spaced bins) of durations.

    Parameters
    ----------
    durations : array-like
    nbins : int
    log_bins : bool

    Returns
    -------
    dict with 'counts', 'edges', 'bin_centers', 'pdf'
    """
    d = np.array(durations, float)

    if log_bins:
        d = d[d > 0]
        if len(d) == 0:
            return {"counts": [], "edges": [], "bin_centers": [], "pdf": []}
        edges = np.logspace(np.log10(d.min()), np.log10(d.max()), nbins + 1)
        bin_centers = np.sqrt(edges[:-1] * edges[1:])
    else:
        edges = np.linspace(d.min(), d.max(), nbins + 1)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])

    counts, edges = np.histogram(d, bins=edges)
    bin_widths = np.diff(edges)
    pdf = counts / (bin_widths * len(d))

    return {
        "counts": counts,
        "edges": edges,
        "bin_centers": bin_centers,
        "pdf": pdf,
    }


def bootstrap_ci(data, statistic="mean", n_boot=5000):
    """
    Bootstrap confidence interval.

    Parameters
    ----------
    data : array-like
    statistic : {"mean", "median"}
    n_boot : int

    Returns
    -------
    ci : ndarray, shape (2,)
        2.5 and 97.5 percentiles.
    """
    data = np.asarray(data)
    stats = []
    for _ in range(n_boot):
        boot = np.random.choice(data, size=len(data), replace=True)
        if statistic == "mean":
            stats.append(np.mean(boot))
        elif statistic == "median":
            stats.append(np.median(boot))
        else:
            raise ValueError("Unknown statistic: %s" % statistic)

    return np.percentile(stats, [2.5, 97.5])
