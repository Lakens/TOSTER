import numpy as np
import pandas as pd
import warnings
from .rbs import rbs_calc, rb_to_cstat, rb_to_odds
from .boot_ci import basic_ci, perc_ci, stud_ci

def boot_ses_calc(x, y=None, paired=False, ses="rb", alpha=0.05, boot_ci="basic", R=1999, mu=0):

    x = np.array(x)
    if y is not None:
        y = np.array(y)

    raw_ses_val = rbs_calc(x, y, mu, paired)

    boots = np.zeros(R)
    if paired:
        data = x if y is None else x - y
        data = data[~np.isnan(data)]
        for i in range(R):
            sampler = np.random.choice(data, len(data), replace=True)
            boots[i] = rbs_calc(sampler, mu=mu, paired=True)
    else:
        if y is None:
            data = x[~np.isnan(x)]
            for i in range(R):
                sampler = np.random.choice(data, len(data), replace=True)
                boots[i] = rbs_calc(sampler, mu=mu, paired=True)
        else:
            x_clean = x[~np.isnan(x)]
            y_clean = y[~np.isnan(y)]
            for i in range(R):
                x_boot = np.random.choice(x_clean, len(x_clean), replace=True)
                y_boot = np.random.choice(y_clean, len(y_clean), replace=True)
                boots[i] = rbs_calc(x_boot, y_boot, mu=mu, paired=False)

    z_boots = np.arctanh(boots)
    z_raw = np.arctanh(raw_ses_val)

    if boot_ci == "perc":
        zci = perc_ci(z_boots, alpha)
    elif boot_ci == "basic":
        zci = basic_ci(z_boots, z_raw, alpha)
    elif boot_ci == "stud":
        warnings.warn("Studentized CI for rank-biserial correlation is not fully implemented, using percentile CI instead.")
        zci = perc_ci(z_boots, alpha)

    rci = np.tanh(zci)

    if ses == 'rb':
        estimate = raw_ses_val
        ci = rci
        boots2 = boots
    elif ses == 'cstat':
        estimate = rb_to_cstat(raw_ses_val)
        ci = rb_to_cstat(rci)
        boots2 = rb_to_cstat(boots)
    elif ses == 'odds':
        estimate = rb_to_odds(raw_ses_val)
        ci = rb_to_odds(rci)
        boots2 = rb_to_odds(boots)
    elif ses == 'logodds':
        estimate = np.log(rb_to_odds(raw_ses_val))
        ci = np.log(rb_to_odds(rci))
        boots2 = np.log(rb_to_odds(boots))

    effsize = pd.DataFrame({
        'estimate': estimate,
        'bias': estimate - np.median(boots2),
        'SE': np.std(boots2),
        'lower.ci': ci[0],
        'upper.ci': ci[1],
        'conf.level': 1-alpha,
        'boot_ci': boot_ci
    }, index=[ses])

    return effsize
