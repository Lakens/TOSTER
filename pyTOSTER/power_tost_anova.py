import numpy as np
from scipy import stats
from scipy.optimize import brentq
import warnings

def power_tost_anova(alpha=0.05, df1=None, df2=None, eqbound=None, power=None):

    if sum(x is None for x in [alpha, df1, df2, eqbound, power]) != 1:
        raise ValueError("Exactly one of alpha, df1, df2, eqbound, or power must be None")

    def calc_power(alpha, df1, df2, eqbound):
        if alpha is None or df1 is None or df2 is None or eqbound is None:
            return np.nan
        f2 = eqbound / (1 - eqbound)
        lambda_ = f2 * (df1 + df2 + 1)
        F_crit = stats.f.ppf(alpha, df1, df2, nc=lambda_)
        return stats.f.cdf(F_crit, df1, df2)

    if power is None:
        power = calc_power(alpha, df1, df2, eqbound)
    elif df2 is None:
        func = lambda df2: calc_power(alpha, df1, df2, eqbound) - power
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df2 = brentq(func, a=4, b=1e5)
    elif eqbound is None:
        func = lambda eqbound: calc_power(alpha, df1, df2, eqbound) - power
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            eqbound = brentq(func, a=0.001, b=0.999)
    elif alpha is None:
        func = lambda alpha: calc_power(alpha, df1, df2, eqbound) - power
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            alpha = brentq(func, a=1e-10, b=1-1e-10)
    elif df1 is None:
        func = lambda df1: calc_power(alpha, df1, df2, eqbound) - power
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df1 = brentq(func, a=1, b=100)

    return {
        'alpha': alpha,
        'df1': df1,
        'df2': df2,
        'eqbound': eqbound,
        'power': power,
        'method': "Power Analysis for F-test Equivalence Testing"
    }
