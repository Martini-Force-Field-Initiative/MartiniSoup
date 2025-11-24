import numpy as np
from lmfit import Model

def _exp_survival(t, koff):
    return np.exp(-koff * t)

class SingleExponentialModel:
    """
    Fit S(t) = exp(-koff * t)
    """

    def __init__(self, times, survival):
        self.times = np.asarray(times, dtype=float)
        self.survival = np.asarray(survival, dtype=float)

    def fit(self, koff_init=1.0):
        """
        Returns
        -------
        lmfit.ModelResult
        """
        model = Model(_exp_survival)
        params = model.make_params(koff=koff_init)
        result = model.fit(self.survival, params, t=self.times)
        return result
