import numpy as np
from scipy import stats

def hedge_j(df):
    if df < 2: return np.nan
    return np.exp(stats.lgamma(df / 2) - np.log(np.sqrt(df / 2)) - stats.lgamma((df - 1) / 2))

def d_est_one(n, m, sd, smd_type='d'):
    if sd == 0: return {'d': np.inf, 'se': np.nan}
    d = m / sd
    df = n - 1
    j = hedge_j(df) if smd_type == 'g' else 1
    d = d * j
    se = np.sqrt(1/n + (d**2 / (2*n))) if n > 0 else np.nan
    return {'d': d, 'se': se}

def d_est_pair(n, m1, m2, sd1, sd2, r12, smd_type='d', rm_correction=False):
    sd_diff = np.sqrt(sd1**2 + sd2**2 - 2 * r12 * sd1 * sd2)
    if sd_diff == 0: return {'d': np.inf, 'se': np.nan}
    if rm_correction:
        denom = sd_diff / np.sqrt(2 * (1 - r12)) if (1-r12) > 0 else 0
    else:
        denom = sd_diff
    if denom == 0: return {'d': np.inf, 'se': np.nan}
    d = (m1 - m2) / denom
    df = n - 1
    j = hedge_j(df) if smd_type == 'g' else 1
    d = d * j
    se = np.sqrt(1/n + (d**2 / (2*n))) if n > 0 else np.nan
    return {'d': d, 'se': se}

def d_est_ind(n1, n2, m1, m2, sd1, sd2, smd_type='d', var_equal=False, glass=None):
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

    j = hedge_j(df) if smd_type == 'g' else 1
    d = d * j

    se = np.sqrt((n1+n2)/(n1*n2) + (d**2 / (2*(n1+n2)))) if n1 > 0 and n2 > 0 else np.nan

    return {'d': d, 'se': se}

def pool_sd(x, y):
    n1, n2 = len(x), len(y)
    if n1 + n2 <= 2: return np.nan
    sd1, sd2 = np.std(x, ddof=1), np.std(y, ddof=1)
    return np.sqrt((((n1 - 1) * (sd1**2)) + ((n2 - 1) * (sd2**2))) / ((n1 + n2) - 2))
