

"""
LMfit-based kinetic models for survival-curve fitting.

Models included:
    • Single exponential S(t) = exp(-koff * t)
    • Weibull survival
    • Double exponential mixture
"""

import numpy as np
from lmfit import Model

# ------------------------------------------------------------------
# MODELS

def exp_survival(t, koff):
    return np.exp(-koff * t)


def weibull_survival(t, k, beta):
    return np.exp(-(k * t) ** beta)


def double_exp_survival(t, A, k1, k2):
    return A * np.exp(-k1 * t) + (1 - A) * np.exp(-k2 * t)


# Create LMfit model objects
exp_model = Model(exp_survival)
weibull_mdl = Model(weibull_survival)
double_mdl = Model(double_exp_survival)


class KineticModels:
    """
    Wrapper class for fitting kinetic survival functions using LMfit.
    """

    @staticmethod
    def fit_exponential(time, survival):
        """Fit single-exponential survival curve."""
        return exp_model.fit(survival, t=time, koff=1e-3)

    @staticmethod
    def fit_weibull(time, survival):
        """Fit Weibull survival curve."""
        return weibull_mdl.fit(survival, t=time, k=1e-3, beta=1.0)

    @staticmethod
    def fit_double_exponential(time, survival):
        """Fit a two-state double-exponential survival model."""
        return double_mdl.fit(
            survival,
            t=time,
            A=0.5,
            k1=1e-2,
            k2=1e-4,
        )
