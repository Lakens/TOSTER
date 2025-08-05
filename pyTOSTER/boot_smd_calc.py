import numpy as np
from scipy import stats
import pandas as pd
import warnings

def _hedge_j(df):
    if df < 2: return np.nan
    return np.exp(stats.lgamma(df / 2) - np.log(np.sqrt(df / 2)) - stats.lgamma((df - 1) / 2))

def _d_est_one(n, m, sd, smd_type='d'):
    if sd == 0: return {'d': np.inf, 'se': np.nan}
    d = m / sd
    df = n - 1
    j = _hedge_j(df) if smd_type == 'g' else 1
    d = d * j
    se = np.sqrt(1/n + (d**2 / (2*n))) if n > 0 else np.nan
    return {'d': d, 'se': se}

def _d_est_pair(n, m1, m2, sd1, sd2, r12, smd_type='d', rm_correction=False):
    sd_diff = np.sqrt(sd1**2 + sd2**2 - 2 * r12 * sd1 * sd2)
    if sd_diff == 0: return {'d': np.inf, 'se': np.nan}
    if rm_correction:
        denom = sd_diff / np.sqrt(2 * (1 - r12)) if (1-r12) > 0 else 0
    else:
        denom = sd_diff
    if denom == 0: return {'d': np.inf, 'se': np.nan}
    d = (m1 - m2) / denom
    df = n - 1
    j = _hedge_j(df) if smd_type == 'g' else 1
    d = d * j
    se = np.sqrt(1/n + (d**2 / (2*n))) if n > 0 else np.nan
    return {'d': d, 'se': se}

def _d_est_ind(n1, n2, m1, m2, sd1, sd2, smd_type='d', var_equal=False, glass=None):
    if glass == 'glass1':
        sd_denom = sd1
    elif glass == 'glass2':
        sd_denom = sd2
    else:
        if var_equal:
            sd_denom = np.sqrt((((n1 - 1) * sd1**2) + ((n2 - 1) * sd2**2)) / (n1 + n2 - 2)) if n1 + n2 > 2 else 0
        else:
            sd_denom = np.sqrt((sd1**2 + sd2**2) / 2)

    if sd_denom == 0: return {'d': np.inf, 'se': np.nan}
    d = (m1 - m2) / sd_denom

    if var_equal:
        df = n1 + n2 - 2
    else:
        df = (sd1**2/n1 + sd2**2/n2)**2 / (((sd1**2/n1)**2/(n1-1)) + ((sd2**2/n2)**2/(n2-1))) if n1 > 1 and n2 > 1 else 0

    j = _hedge_j(df) if smd_type == 'g' else 1
    d = d * j

    se = np.sqrt((n1+n2)/(n1*n2) + (d**2 / (2*(n1+n2)))) if n1 > 0 and n2 > 0 else np.nan

    return {'d': d, 'se': se}

def _smd_calc(x, y=None, paired=False, var_equal=False, bias_correction=True, rm_correction=False, glass=None, mu=0):
    smd_type = 'g' if bias_correction else 'd'
    x = np.array(x)
    if y is not None: y = np.array(y)

    if y is None: # one-sample
        x = x[~np.isnan(x)]
        res = _d_est_one(len(x), np.mean(x) - mu, np.std(x, ddof=1), smd_type)
    elif paired:
        mask = ~np.isnan(x) & ~np.isnan(y)
        x, y = x[mask], y[mask]
        if len(x) < 2: return pd.DataFrame({'estimate': np.nan, 'SE': np.nan}, index=[smd_type])
        r12, _ = stats.pearsonr(x, y)
        res = _d_est_pair(len(x), np.mean(x), np.mean(y), np.std(x, ddof=1), np.std(y, ddof=1), r12, smd_type, rm_correction)
    else: # independent
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]
        res = _d_est_ind(len(x), len(y), np.mean(x), np.mean(y), np.std(x, ddof=1), np.std(y, ddof=1), smd_type, var_equal, glass)

    return pd.DataFrame({'estimate': res['d'], 'SE': res['se']}, index=[smd_type])

def boot_smd_calc(x, y=None, paired=False, var_equal=False, alpha=0.05, bias_correction=True, rm_correction=False, glass=None, boot_ci="stud", R=1999, mu=0):

    raw_smd = _smd_calc(x, y, paired, var_equal, bias_correction, rm_correction, glass, mu)

    boots = np.zeros(R)
    boots_se = np.zeros(R)

    if y is None:
        x = np.array(x)[~np.isnan(np.array(x))]
        for i in range(R):
            sampler = np.random.choice(x, len(x), replace=True)
            res_boot = _smd_calc(sampler, None, False, False, bias_correction, False, None, mu)
            boots[i] = res_boot['estimate'].iloc[0]
            boots_se[i] = res_boot['SE'].iloc[0]
    elif paired:
        mask = ~np.isnan(np.array(x)) & ~np.isnan(np.array(y))
        x, y = np.array(x)[mask], np.array(y)[mask]
        for i in range(R):
            sampler = np.random.choice(len(x), len(x), replace=True)
            res_boot = _smd_calc(x[sampler], y[sampler], True, False, bias_correction, rm_correction, None, mu)
            boots[i] = res_boot['estimate'].iloc[0]
            boots_se[i] = res_boot['SE'].iloc[0]
    else:
        x = np.array(x)[~np.isnan(np.array(x))]
        y = np.array(y)[~np.isnan(np.array(y))]
        for i in range(R):
            x_boot = np.random.choice(x, len(x), replace=True)
            y_boot = np.random.choice(y, len(y), replace=True)
            res_boot = _smd_calc(x_boot, y_boot, False, var_equal, bias_correction, False, glass, mu)
            boots[i] = res_boot['estimate'].iloc[0]
            boots_se[i] = res_boot['SE'].iloc[0]

    if boot_ci == "stud":
        ci = _stud_ci(boots, boots_se, raw_smd['SE'].iloc[0], raw_smd['estimate'].iloc[0], alpha)
    elif boot_ci == "basic":
        ci = _basic_ci(boots, raw_smd['estimate'].iloc[0], alpha)
    else: # perc
        ci = _perc_ci(boots, alpha)

    effsize = pd.DataFrame({
        'estimate': raw_smd['estimate'].iloc[0],
        'bias': raw_smd['estimate'].iloc[0] - np.nanmedian(boots),
        'SE': np.nanstd(boots),
        'lower.ci': ci[0],
        'upper.ci': ci[1],
        'conf.level': 1-alpha,
        'boot_ci': boot_ci
    }, index=[raw_smd.index[0]])

    return effsize

def _perc_ci(boots_est, alpha):
    return np.quantile(boots_est, [alpha/2, 1-alpha/2])

def _basic_ci(boots_est, t0, alpha):
    conf = 1 - alpha
    qq = np.quantile(boots_est, [(1 - conf) / 2, (1 + conf) / 2])
    return 2 * t0 - qq[::-1]

def _stud_ci(boots_est, boots_se, se0, t0, alpha):
    conf = 1 - alpha
    z = (boots_est - t0) / boots_se
    z = z[np.isfinite(z)]
    qq = np.quantile(z, [(1 - conf) / 2, (1 + conf) / 2])
    return t0 - se0 * qq[::-1]
