from .survival import SurvivalAnalysis
from .histogram import HistogramAnalysis
from .models import SingleExponentialModel

class ResidenceAnalysis:
    """
    High-level interface for performing analyses on the
    type_agg_named dictionary.
    """

    def __init__(self, type_agg_named):
        """
        Parameters
        ----------
        type_agg_named : dict[str, list[float]]
            Mapping molecule type -> list of residence durations
        """
        self.data = type_agg_named

    def types(self):
        """Return available metabolite types."""
        return list(self.data.keys())

    def durations(self, typename):
        """Return durations for a given metabolite type."""
        if typename is not None:
            return self.data[typename]
        else:
            raise KeyError('Must specify molecule')

    # --------------- Survival Curves ----------------

    def survival(self, typename):
        sa = SurvivalAnalysis(self.durations(typename))
        return sa.compute_survival()

    def survival_CI(self, typename, n_boot=200):
        sa = SurvivalAnalysis(self.durations(typename))
        sa.compute_survival()
        return sa.bootstrap_CI(n_boot=n_boot)

    # --------------- Histogram ----------------

    def histogram(self, typename, nbins=30, log_bins=False):
        ha = HistogramAnalysis(self.durations(typename))
        return ha.compute(nbins=nbins, log_bins=log_bins)

    # --------------- Kinetic model fitting ----------------

    def fit_exponential(self, typename, koff_init=1.0):
        sa = SurvivalAnalysis(self.durations(typename))
        t, S = sa.compute_survival()
        model = SingleExponentialModel(t, S)
        return model.fit(koff_init=koff_init)
