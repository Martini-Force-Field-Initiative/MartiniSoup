
"""
Analysis utilities for residence-time data.

Includes:
    • Histogram analysis (linear or log-binned)
    • Survival curve computation
    • Hazard function estimation
    • Bootstrap confidence intervals
"""

import numpy as np
from collections import Counter
from .utils import binned_histogram, bootstrap_ci


class ResidenceAnalysis:
    """
    Perform analysis operations on aggregated residence-time data.

    Parameters
    ----------
    residences : MetaboliteResidences
        Aggregated residence-time container.
    """

    def __init__(self, residences):
        self.data = residences

    # ------------------------------------------------------------------
    def histogram(self, moltype_id, log_bins=True, nbins=20):
        """
        Compute a histogram of residence durations.

        Parameters
        ----------
        moltype_id : int
            Molecule type ID.
        log_bins : bool, default True
            Use logarithmic binning.
        nbins : int
            Number of bins.

        Returns
        -------
        dict
            Contains 'counts', 'edges', 'bin_centers', and 'pdf'.
        """
        durations = self.data.type_data[moltype_id]
        return binned_histogram(durations, nbins=nbins, log_bins=log_bins)

    # ------------------------------------------------------------------
    def survival_curve(self, moltype_id):
        """
        Construct the empirical survival function S(t).

        Parameters
        ----------
        moltype_id : int

        Returns
        -------
        dict with keys:
            'time' : array
            'survival' : array
        """
        durations = np.array(self.data.type_data[moltype_id], float)
        event_count = Counter(durations)
        times = np.array(sorted(event_count.keys()))

        S = []
        n_total = len(durations)
        alive = n_total

        for t in times:
            S.append(alive / n_total)
            alive -= event_count[t]

        return {"time": times, "survival": np.array(S)}

    # ------------------------------------------------------------------
    def hazard(self, survival_curve):
        """
        Estimate the hazard function h(t) from the survival curve.

        Parameters
        ----------
        survival_curve : dict
            Output of `survival_curve()`.

        Returns
        -------
        dict with keys 'time' and 'hazard'.
        """
        t = survival_curve["time"]
        S = survival_curve["survival"]
        hazard = -np.gradient(np.log(S), t)
        return {"time": t, "hazard": hazard}

    # ------------------------------------------------------------------
    def confidence_interval(self, moltype_id, statistic="mean", n_boot=5000):
        """
        Bootstrap confidence interval for a statistic of residence durations.

        Parameters
        ----------
        moltype_id : int
        statistic : {"mean", "median"}
            Statistic to compute.
        n_boot : int
            Number of bootstrap resamples.

        Returns
        -------
        ndarray
            2.5 and 97.5 percentile confidence interval.
        """
        durations = np.array(self.data.type_data[moltype_id])
        return bootstrap_ci(durations, statistic=statistic, n_boot=n_boot)
