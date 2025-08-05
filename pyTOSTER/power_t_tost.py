import numpy as np
from scipy.optimize import brentq
from .power_tost_functions import pow_t_tost
import warnings

def power_t_tost(n=None, delta=0, sd=1, eqb=None, low_eqbound=None, high_eqbound=None,
                 alpha=None, power=None, type="two.sample"):

    if eqb is not None:
        if isinstance(eqb, (int, float)):
            high_eqbound = abs(eqb)
            low_eqbound = -abs(eqb)
        else:
            high_eqbound = max(eqb)
            low_eqbound = min(eqb)

    if low_eqbound is None or high_eqbound is None:
        raise ValueError("Equivalence bounds must be provided.")

    if low_eqbound > delta or high_eqbound < delta:
        raise ValueError("True mean difference greater than bounds. TOST power calculation not possible.")

    if sum(x is None for x in [n, power, alpha]) != 1:
        raise ValueError("Exactly one of n, power, or alpha must be NULL")

    if n is None:
        func = lambda n: pow_t_tost(alpha=alpha, theta1=low_eqbound, theta2=high_eqbound,
                                     theta0=delta, sd=sd, n=n, type=type) - power
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            n = brentq(func, a=2 + 1e-10, b=1e+09)
    elif power is None:
        power = pow_t_tost(alpha=alpha, theta1=low_eqbound, theta2=high_eqbound,
                           theta0=delta, sd=sd, n=n, type=type)
    elif alpha is None:
        func = lambda alpha: pow_t_tost(alpha=alpha, theta1=low_eqbound, theta2=high_eqbound,
                                         theta0=delta, sd=sd, n=n, type=type) - power
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            alpha = brentq(func, a=1e-10, b=1 - 1e-10)

    return {
        'power': power,
        'alpha': alpha,
        'n': n,
        'delta': delta,
        'sd': sd,
        'low_eqbound': low_eqbound,
        'high_eqbound': high_eqbound
    }
