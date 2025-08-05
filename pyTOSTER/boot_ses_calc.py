import numpy as np
from scipy import stats
import pandas as pd
from scipy.special import expit, logit
import warnings

def _ranktransform(x, sign=False, method='average'):
    if np.all(np.isnan(x)):
        return x
    if len(np.unique(x)) == 1:
        return x
    if sign:
        zeros = x == 0
        out = np.full(len(x), np.nan)
        if np.any(~zeros):
            out[~zeros] = np.sign(x[~zeros]) * stats.rankdata(np.abs(x[~zeros]), method=method)
        return out
    else:
        return stats.rankdata(x, method=method)

def _rbs_calc(x, y=None, mu=0, paired=False):
    if paired:
        if y is None:
            z = x - mu
        else:
            z = (x - y) - mu
        z = z[~np.isnan(z)]
        if len(z) == 0: return 0.0
        abs_z = np.abs(z)
        RR = -1 * _ranktransform(abs_z) * np.sign(z)
        Rplus = np.sum(RR[RR > 0])
        Rminus = np.sum(np.abs(RR[RR < 0]))
        Tee = min(Rplus, Rminus)
        n = len(RR)
        if n == 0: return 0.0
        if Rplus >= Rminus:
            rho = -4 * abs((Tee - (Rplus + Rminus) / 2) / n / (n + 1))
        else:
            rho = 4 * abs((Tee - (Rplus + Rminus) / 2) / n / (n + 1))
        return rho
    else:
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]
        if len(x) == 0 or len(y) == 0: return 0.0
        Ry = _ranktransform(np.concatenate([x - mu, y]))
        n1 = len(x)
        n2 = len(y)
        S = n1 * n2
        U1 = np.sum(Ry[:n1]) - n1 * (n1 + 1) / 2
        # Cliff's delta, equivalent to rank-biserial correlation for independent samples
        return (U1 / S) * 2 - 1

def _rb_to_cstat(x):
    return (x + 1) / 2

def _rb_to_odds(x):
    return np.exp(logit(_rb_to_cstat(x)))

def boot_ses_calc(x, y=None, paired=False, ses="rb", alpha=0.05, boot_ci="basic", R=1999, mu=0):

    x = np.array(x)
    if y is not None:
        y = np.array(y)

    raw_ses_val = _rbs_calc(x, y, mu, paired)

    boots = np.zeros(R)
    if paired:
        data = x if y is None else x - y
        data = data[~np.isnan(data)]
        for i in range(R):
            sampler = np.random.choice(data, len(data), replace=True)
            boots[i] = _rbs_calc(sampler, mu=mu, paired=True)
    else:
        if y is None:
            data = x[~np.isnan(x)]
            for i in range(R):
                sampler = np.random.choice(data, len(data), replace=True)
                boots[i] = _rbs_calc(sampler, mu=mu, paired=True)
        else:
            x_clean = x[~np.isnan(x)]
            y_clean = y[~np.isnan(y)]
            for i in range(R):
                x_boot = np.random.choice(x_clean, len(x_clean), replace=True)
                y_boot = np.random.choice(y_clean, len(y_clean), replace=True)
                boots[i] = _rbs_calc(x_boot, y_boot, mu=mu, paired=False)

    z_boots = np.arctanh(boots)
    z_raw = np.arctanh(raw_ses_val)

    if boot_ci == "perc":
        zci = np.quantile(z_boots, [alpha/2, 1-alpha/2])
    elif boot_ci == "basic":
        zci = 2 * z_raw - np.quantile(z_boots, [1-alpha/2, alpha/2])
    elif boot_ci == "stud":
        warnings.warn("Studentized CI for rank-biserial correlation is not fully implemented, using percentile CI instead.")
        zci = np.quantile(z_boots, [alpha/2, 1-alpha/2])

    rci = np.tanh(zci)

    if ses == 'rb':
        estimate = raw_ses_val
        ci = rci
        boots2 = boots
    elif ses == 'cstat':
        estimate = _rb_to_cstat(raw_ses_val)
        ci = _rb_to_cstat(rci)
        boots2 = _rb_to_cstat(boots)
    elif ses == 'odds':
        estimate = _rb_to_odds(raw_ses_val)
        ci = _rb_to_odds(rci)
        boots2 = _rb_to_odds(boots)
    elif ses == 'logodds':
        estimate = np.log(_rb_to_odds(raw_ses_val))
        ci = np.log(_rb_to_odds(rci))
        boots2 = np.log(_rb_to_odds(boots))

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
