import numpy as np
from scipy import stats

def power_tost_two(alpha=None, statistical_power=None, N=None, low_eqbound_d=None, high_eqbound_d=None):
    if sum(x is None for x in [N, statistical_power, low_eqbound_d, high_eqbound_d]) != 1:
        # special case for both eqbounds being None
        if not (low_eqbound_d is None and high_eqbound_d is None and N is not None and statistical_power is not None):
            raise ValueError("Exactly one of N, statistical_power, or eqbounds must be None")

    if N is None:
        nt1 = 2 * (stats.norm.ppf(1 - alpha) + stats.norm.ppf(1 - (1 - statistical_power) / 2))**2 / low_eqbound_d**2
        nt2 = 2 * (stats.norm.ppf(1 - alpha) + stats.norm.ppf(1 - (1 - statistical_power) / 2))**2 / high_eqbound_d**2
        N = max(nt1, nt2)
        return np.ceil(N)
    elif statistical_power is None:
        power1 = 2 * (stats.norm.cdf(abs(low_eqbound_d) * np.sqrt(N / 2) - stats.norm.ppf(1 - alpha)) + stats.norm.cdf(-abs(low_eqbound_d) * np.sqrt(N / 2) - stats.norm.ppf(1 - alpha))) - 1
        power2 = 2 * (stats.norm.cdf(abs(high_eqbound_d) * np.sqrt(N / 2) - stats.norm.ppf(1 - alpha)) + stats.norm.cdf(-abs(high_eqbound_d) * np.sqrt(N / 2) - stats.norm.ppf(1 - alpha))) - 1
        statistical_power = min(power1, power2)
        return max(0, statistical_power)
    elif low_eqbound_d is None and high_eqbound_d is None:
        eqbound = np.sqrt(2 * (stats.norm.ppf(1 - alpha) + stats.norm.ppf(1 - (1 - statistical_power) / 2))**2 / N)
        return -eqbound, eqbound
    else:
        raise ValueError("Invalid combination of arguments.")

def power_tost_two_raw(alpha=None, statistical_power=None, N=None, sdpooled=None, low_eqbound=None, high_eqbound=None):
    if sdpooled is None:
        raise ValueError("sdpooled must be provided for raw power calculation")

    low_eqbound_d = None
    high_eqbound_d = None

    if low_eqbound is not None:
        low_eqbound_d = low_eqbound / sdpooled
    if high_eqbound is not None:
        high_eqbound_d = high_eqbound / sdpooled

    res = power_tost_two(alpha, statistical_power, N, low_eqbound_d, high_eqbound_d)

    if isinstance(res, tuple):
        return res[0] * sdpooled, res[1] * sdpooled
    else:
        return res
