import numpy as np
from statsmodels.stats.power import tt_ind_solve_power, tt_solve_power

def pow_t_tost(alpha=0.05, theta1=None, theta2=None, theta0=0, sd=1, n=None, type="two.sample"):

    if type == "two.sample":
        if n is None:
            raise ValueError("n must be provided for two-sample power analysis")
        # Simplified: assumes equal n
        power1 = tt_ind_solve_power(effect_size=(theta0 - theta1) / sd, nobs1=n, alpha=alpha, power=None, ratio=1.0, alternative='larger')
        power2 = tt_ind_solve_power(effect_size=(theta2 - theta0) / sd, nobs1=n, alpha=alpha, power=None, ratio=1.0, alternative='larger')
    else: # one-sample or paired
        if n is None:
            raise ValueError("n must be provided for one-sample/paired power analysis")
        power1 = tt_solve_power(effect_size=(theta0 - theta1) / sd, nobs=n, alpha=alpha, power=None, alternative='larger')
        power2 = tt_solve_power(effect_size=(theta2 - theta0) / sd, nobs=n, alpha=alpha, power=None, alternative='larger')

    return power1 + power2 - 1
