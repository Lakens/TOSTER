import numpy as np
from scipy import stats
import warnings

def se_dz(smd, n):
    if n <= 1: return np.nan
    return np.sqrt(1/n + (smd**2 / (2*n)))

def se_ds(smd, n):
    if isinstance(n, int):
        n = [n, n]
    if n[0] + n[1] <= 2: return np.nan
    return np.sqrt((n[0] + n[1]) / (n[0] * n[1]) + smd**2 / (2 * (n[0] + n[1])))

def pool_sd(x, y):
    n1, n2 = len(x), len(y)
    if n1 + n2 <= 2: return np.nan
    sd1, sd2 = np.std(x, ddof=1), np.std(y, ddof=1)
    return np.sqrt((((n1 - 1) * (sd1**2)) + ((n2 - 1) * (sd2**2))) / ((n1 + n2) - 2))


def boot_compare_smd(x1, y1=None, x2=None, y2=None,
                     null=0, paired=False,
                     alternative='two.sided',
                     R=1999, alpha=0.05):
    """
    Comparing Standardized Mean Differences (SMDs) Between Independent Studies with Bootstrapping.
    """
    if not isinstance(x1, (list, np.ndarray)):
        raise TypeError("x1 must be a list or numpy array.")

    # Data preparation
    if paired:
        if y1 is None and y2 is None:
            df1 = np.array(x1)
            df2 = np.array(x2)
        else:
            df1 = np.array(x1) - np.array(y1)
            df2 = np.array(x2) - np.array(y2)
        df1 = df1[~np.isnan(df1)]
        df2 = df2[~np.isnan(df2)]
    else:
        if y1 is None or y2 is None:
            raise ValueError("For independent samples, y1 and y2 must be provided.")
        x1, y1 = np.array(x1), np.array(y1)
        x2, y2 = np.array(x2), np.array(y2)
        x1 = x1[~np.isnan(x1)]
        y1 = y1[~np.isnan(y1)]
        x2 = x2[~np.isnan(x2)]
        y2 = y2[~np.isnan(y2)]

    # Bootstrap loop
    smd1_vec = np.zeros(R)
    smd2_vec = np.zeros(R)

    for i in range(R):
        if paired:
            df1_boot = np.random.choice(df1, len(df1), replace=True)
            df2_boot = np.random.choice(df2, len(df2), replace=True)
            smd1_vec[i] = np.mean(df1_boot) / np.std(df1_boot, ddof=1)
            smd2_vec[i] = np.mean(df2_boot) / np.std(df2_boot, ddof=1)
        else:
            x1_boot = np.random.choice(x1, len(x1), replace=True)
            y1_boot = np.random.choice(y1, len(y1), replace=True)
            x2_boot = np.random.choice(x2, len(x2), replace=True)
            y2_boot = np.random.choice(y2, len(y2), replace=True)

            md1_boot = np.mean(x1_boot) - np.mean(y1_boot)
            sd1_boot = pool_sd(x1_boot, y1_boot)
            smd1_vec[i] = md1_boot / sd1_boot if sd1_boot != 0 else 0

            md2_boot = np.mean(x2_boot) - np.mean(y2_boot)
            sd2_boot = pool_sd(x2_boot, y2_boot)
            smd2_vec[i] = md2_boot / sd2_boot if sd2_boot != 0 else 0

    d_diff_vec = smd1_vec - smd2_vec

    # P-value and CI
    if alternative in ['equivalence', 'minimal.effect']:
        if isinstance(null, (int, float)):
            low_eqbound, high_eqbound = -abs(null), abs(null)
        else:
            low_eqbound, high_eqbound = min(null), max(null)

        if alternative == 'equivalence':
            p1 = np.sum(d_diff_vec >= high_eqbound) / R
            p2 = np.sum(d_diff_vec <= low_eqbound) / R
            p_value = max(p1, p2)
        else: # minimal.effect
            p1 = np.sum(d_diff_vec <= high_eqbound) / R
            p2 = np.sum(d_diff_vec >= low_eqbound) / R
            p_value = min(p1, p2)

        ci_level = 1 - 2 * alpha
    else:
        if alternative == 'two.sided':
            phat = (np.sum(d_diff_vec < null) + 0.5 * np.sum(d_diff_vec == null)) / R
            p_value = 2 * min(phat, 1 - phat)
        elif alternative == 'greater':
            p_value = np.sum(d_diff_vec <= null) / R
        elif alternative == 'less':
            p_value = np.sum(d_diff_vec >= null) / R

        ci_level = 1 - alpha

    conf_int = np.quantile(d_diff_vec, [(1 - ci_level) / 2, 1 - (1 - ci_level) / 2])

    # Initial estimates
    if paired:
        smd1 = np.mean(df1) / np.std(df1, ddof=1)
        smd2 = np.mean(df2) / np.std(df2, ddof=1)
    else:
        md1 = np.mean(x1) - np.mean(y1)
        sd1 = pool_sd(x1, y1)
        smd1 = md1 / sd1 if sd1 != 0 else 0
        md2 = np.mean(x2) - np.mean(y2)
        sd2 = pool_sd(x2, y2)
        smd2 = md2 / sd2 if sd2 != 0 else 0

    estimate = smd1 - smd2

    # Standard error of the difference
    if paired:
        se1 = se_dz(smd1, len(df1))
        se2 = se_dz(smd2, len(df2))
    else:
        se1 = se_ds(smd1, [len(x1), len(y1)])
        se2 = se_ds(smd2, [len(x2), len(y2)])

    stderr = np.sqrt(se1**2 + se2**2) if not np.isnan(se1) and not np.isnan(se2) else np.nan
    z_stat = (estimate - null) / stderr if not isinstance(null, list) and stderr != 0 else np.nan

    results = {
        'statistic': z_stat,
        'p.value': p_value,
        'conf.int': conf_int,
        'estimate': estimate,
        'null.value': null,
        'alternative': alternative,
        'method': 'Bootstrapped Differences in SMDs',
        'boot_res': {'smd1': smd1_vec, 'smd2': smd2_vec, 'd_diff': d_diff_vec}
    }

    return results
