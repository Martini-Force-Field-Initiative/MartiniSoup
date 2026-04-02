import numpy as np
from scipy.stats import rankdata

class SurvivalAnalysis:
    """
    Compute empirical survival curves and bootstrap confidence intervals.
    """

    def __init__(self, durations):
        """
        Parameters
        ----------
        durations : array-like
            List of residence times for a single molecule type.
        """
        self.durations = np.asarray(durations, dtype=float)
        self.times = None
        self.survival = None

    def compute_survival(self):
        """
        Compute the empirical survival function S(t).

        Returns
        -------
        times : ndarray
        S : ndarray
        """
        if len(self.durations) == 0:
            return np.array([]), np.array([])

        # Sort durations
        sorted_T = np.sort(self.durations)
        n = len(sorted_T)

        # Survival = (n - rank + 1) / n
        ranks = rankdata(sorted_T, method="ordinal")
        S = (n - ranks + 1) / n

        self.times = sorted_T
        self.survival = S

        return self.times, self.survival

    def bootstrap_CI(self, n_boot=200):
        """
        Compute bootstrap confidence intervals for survival curves.

        Returns
        -------
        (lower, upper) : tuple of arrays
        """
        if self.times is None or self.survival is None:
            self.compute_survival()

        t = self.times
        n = len(self.durations)
        curves = []

        for _ in range(n_boot):
            sample = np.random.choice(self.durations, size=n, replace=True)
            sa = SurvivalAnalysis(sample)
            tt, ss = sa.compute_survival()

            ss_interp = np.interp(t, tt, ss)
            curves.append(ss_interp)

        curves = np.vstack(curves)
        lower = np.percentile(curves, 2.5, axis=0)
        upper = np.percentile(curves, 97.5, axis=0)

        return lower, upper
