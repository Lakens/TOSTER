import numpy as np
from scipy import stats

def power_tost_paired(alpha=None, statistical_power=None, N=None, low_eqbound_dz=None, high_eqbound_dz=None):
    if sum(x is None for x in [N, statistical_power, low_eqbound_dz, high_eqbound_dz]) != 1:
        # special case for both eqbounds being None
        if not (low_eqbound_dz is None and high_eqbound_dz is None and N is not None and statistical_power is not None):
            raise ValueError("Exactly one of N, statistical_power, or eqbounds must be None")

    if N is None:
        nt1 = (stats.norm.ppf(1 - alpha) + stats.norm.ppf(1 - (1 - statistical_power) / 2))**2 / low_eqbound_dz**2
        nt2 = (stats.norm.ppf(1 - alpha) + stats.norm.ppf(1 - (1 - statistical_power) / 2))**2 / high_eqbound_dz**2
        N = max(nt1, nt2)
        return np.ceil(N)
    elif statistical_power is None:
        power1 = 2 * (stats.norm.cdf(abs(low_eqbound_dz) * np.sqrt(N) - stats.norm.ppf(1 - alpha)) + stats.norm.cdf(-abs(low_eqbound_dz) * np.sqrt(N) - stats.norm.ppf(1 - alpha))) - 1
        power2 = 2 * (stats.norm.cdf(abs(high_eqbound_dz) * np.sqrt(N) - stats.norm.ppf(1 - alpha)) + stats.norm.cdf(-abs(high_eqbound_dz) * np.sqrt(N) - stats.norm.ppf(1 - alpha))) - 1
        statistical_power = min(power1, power2)
        return max(0, statistical_power)
    elif low_eqbound_dz is None and high_eqbound_dz is None:
        eqbound = np.sqrt((stats.norm.ppf(1 - alpha) + stats.norm.ppf(1 - (1 - statistical_power) / 2))**2 / N)
        return -eqbound, eqbound
    else:
        raise ValueError("Invalid combination of arguments.")

def power_tost_paired_raw(alpha=None, statistical_power=None, N=None, sdif=None, low_eqbound=None, high_eqbound=None):
    if sdif is None:
        raise ValueError("sdif must be provided for raw power calculation")

    low_eqbound_dz = None
    high_eqbound_dz = None

    if low_eqbound is not None:
        low_eqbound_dz = low_eqbound / sdif
    if high_eqbound is not None:
        high_eqbound_dz = high_eqbound / sdif

    res = power_tost_paired(alpha, statistical_power, N, low_eqbound_dz, high_eqbound_dz)

    if isinstance(res, tuple):
        return res[0] * sdif, res[1] * sdif
    else:
        return res
