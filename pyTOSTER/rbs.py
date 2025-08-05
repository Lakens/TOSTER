import numpy as np
from scipy import stats
from scipy.special import logit

def ranktransform(x, sign=False, method='average'):
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

def rbs_calc(x, y=None, mu=0, paired=False):
    if paired:
        if y is None:
            z = x - mu
        else:
            z = (x - y) - mu
        z = z[~np.isnan(z)]
        if len(z) == 0: return 0.0
        abs_z = np.abs(z)
        RR = -1 * ranktransform(abs_z) * np.sign(z)
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
        Ry = ranktransform(np.concatenate([x - mu, y]))
        n1 = len(x)
        n2 = len(y)
        S = n1 * n2
        U1 = np.sum(Ry[:n1]) - n1 * (n1 + 1) / 2
        return (U1 / S) * 2 - 1

def rb_to_cstat(x):
    return (x + 1) / 2

def rb_to_odds(x):
    return np.exp(logit(rb_to_cstat(x)))

def rbs(x, y=None, mu=0, conf_level=0.95, paired=False):
    r_rbs = rbs_calc(x, y, mu=mu, paired=paired)

    if paired:
        n = len(x)
        rfSE = np.sqrt((2 * n**3 + 3 * n**2 + n) / 6) / ((n**2 + n) / 2) if n > 0 else np.nan
    else:
        n1 = len(x)
        n2 = len(y)
        rfSE = np.sqrt((n1 + n2 + 1) / (3 * n1 * n2)) if n1 > 0 and n2 > 0 else np.nan

    z = np.arctanh(r_rbs)
    ci = z + stats.norm.ppf([ (1-conf_level)/2, 1-(1-conf_level)/2 ]) * rfSE
    confint = np.tanh(ci)

    return {
        'rbs': r_rbs,
        'conf.int': confint
    }

def np_ses(x, y=None, mu=0, conf_level=0.95, paired=False, ses='rb'):
    rb = rbs(x=x, y=y, mu=mu, conf_level=conf_level, paired=paired)

    if ses == 'rb':
        est = rb['rbs']
        confint = rb['conf.int']
    elif ses == 'cstat':
        est = rb_to_cstat(rb['rbs'])
        confint = rb_to_cstat(rb['conf.int'])
    elif ses == 'odds':
        est = rb_to_odds(rb['rbs'])
        confint = rb_to_odds(rb['conf.int'])
    elif ses == 'logodds':
        est = np.log(rb_to_odds(rb['rbs']))
        confint = np.log(rb_to_odds(rb['conf.int']))

    return {
        'est': est,
        'conf.int': confint
    }
