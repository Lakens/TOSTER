import numpy as np
from scipy import stats
from scipy.optimize import brentq
import warnings

def power_z_cor(n = None,
                rho = None,
                power = None,
                null = 0,
                alpha = None,
                alternative = "two.sided"):

    if alternative == "equivalence":
        if isinstance(null, (int, float)):
            if null == 0:
                raise ValueError("null cannot be zero if alternative is equivalence")
            null = [-abs(null), abs(null)]
        return pow_corr_tost(n=n, rho=rho, power=power, null=null, alpha=alpha)
    else:
        # Simplified version without pwr package
        warnings.warn("Power analysis for standard correlation tests is not fully implemented.")
        return None

def pow_corr_tost(n=None, rho=0, power=None, null=None, alpha=None):
    if sum(x is None for x in [n, power, alpha]) != 1:
        raise ValueError("exactly one of n, power, and alpha must be NULL")

    if n is None:
        za = stats.norm.ppf(1 - alpha)
        zb = stats.norm.ppf(1 - (1 - power) / 2)
        c1 = 0.5 * np.log((1 + min(null)) / (1 - min(null)))
        c2 = 0.5 * np.log((1 + max(null)) / (1 - max(null)))
        nt1 = ((za + zb) / c1)**2 + 3
        nt2 = ((za + zb) / c2)**2 + 3
        return np.ceil(max(nt1, nt2))
    elif power is None:
        se = 1 / np.sqrt(n - 3)
        c1 = 0.5 * np.log((1 + min(null)) / (1 - min(null)))
        c2 = 0.5 * np.log((1 + max(null)) / (1 - max(null)))
        power1 = 2 * (stats.norm.cdf(abs(c1) / se - stats.norm.ppf(1 - alpha)) + stats.norm.cdf(-abs(c1) / se - stats.norm.ppf(1 - alpha))) - 1
        power2 = 2 * (stats.norm.cdf(abs(c2) / se - stats.norm.ppf(1 - alpha)) + stats.norm.cdf(-abs(c2) / se - stats.norm.ppf(1 - alpha))) - 1
        return max(0, min(power1, power2))
    elif alpha is None:
        func = lambda alpha: pow_corr_tost(n=n, rho=rho, power=power, null=null, alpha=alpha) - power
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return brentq(func, a=1e-10, b=1-1e-10)
