import numpy as np
from scipy import stats
from scipy.optimize import brentq
from statsmodels.stats.power import zt_ind_solve_power
import warnings

def power_twoprop(p1=None, p2=None, n=None, null=0, alpha=None, power=None, alternative='two.sided'):

    if alternative == 'equivalence':
        if isinstance(null, (int, float)):
            if null == 0:
                raise ValueError("null cannot be zero for equivalence test")
            null = [-abs(null), abs(null)]
        return pow_prop_tost(p1, p2, n, null, alpha, power)
    else:
        # Simplified version without pwr package
        warnings.warn("Power analysis for standard proportion tests is not fully implemented.")
        return None

def pow_prop_tost(p1, p2, n=None, null=None, alpha=None, power=None):
    if sum(x is None for x in [n, power, alpha]) != 1:
        raise ValueError("exactly one of n, power, and alpha must be NULL")

    if n is None:
        prop_se_sq = (p1 * (1 - p1) + p2 * (1 - p2))
        nt1 = prop_se_sq * ((stats.norm.ppf(1 - alpha) + stats.norm.ppf(1 - (1 - power) / 2)) / (abs(p1 - p2) - abs(min(null))))**2
        nt2 = prop_se_sq * ((stats.norm.ppf(1 - alpha) + stats.norm.ppf(1 - (1 - power) / 2)) / (abs(p1 - p2) - abs(max(null))))**2
        return np.ceil(max(nt1, nt2))
    elif power is None:
        prop_se = np.sqrt((p1 * (1 - p1)) / n + (p2 * (1 - p2)) / n)
        power1 = 2 * (stats.norm.cdf((abs(p1 - p2) - min(null)) / prop_se - stats.norm.ppf(1 - alpha)) + stats.norm.cdf(-(abs(p1 - p2) - min(null)) / prop_se - stats.norm.ppf(1 - alpha))) - 1
        power2 = 2 * (stats.norm.cdf((abs(p1 - p2) - max(null)) / prop_se - stats.norm.ppf(1 - alpha)) + stats.norm.cdf(-(abs(p1 - p2) - max(null)) / prop_se - stats.norm.ppf(1 - alpha))) - 1
        return max(0, min(power1, power2))
    elif alpha is None:
        func = lambda alpha: pow_prop_tost(p1, p2, n, null, alpha, power) - power
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return brentq(func, a=1e-10, b=1-1e-10)
