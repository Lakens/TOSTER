import numpy as np
from scipy import stats
import pandas as pd
import warnings

def _simple_t_test_helper(x, y=None, paired=False, var_equal=False, alternative='two.sided', mu=0, alpha=0.05):
    if y is None:
        res = stats.ttest_1samp(x, popmean=mu)
        stderr = np.std(x, ddof=1) / np.sqrt(len(x))
    elif paired:
        res = stats.ttest_rel(x, y)
        stderr = np.std(x - y, ddof=1) / np.sqrt(len(x))
    else:
        res = stats.ttest_ind(x, y, equal_var=var_equal)
        stderr = np.sqrt(np.var(x, ddof=1)/len(x) + np.var(y, ddof=1)/len(y))

    # This is a simplified htest object
    return {
        'statistic': res.statistic,
        'p.value': res.pvalue,
        'null.value': mu,
        'stderr': stderr
    }

def boot_t_test(x, y=None, var_equal=False, paired=False, alternative='two.sided', mu=0, alpha=0.05, boot_ci="stud", R=1999):

    x = np.array(x)[~np.isnan(np.array(x))]
    if y is not None:
        y = np.array(y)[~np.isnan(np.array(y))]

    null_test = _simple_t_test_helper(x, y, paired, var_equal, alternative, mu, alpha)

    m_vec = np.zeros(R)
    m_se_vec = np.zeros(R)

    if y is None:
        for i in range(R):
            sampler = np.random.choice(x, len(x), replace=True)
            m_vec[i] = np.mean(sampler)
            m_se_vec[i] = np.std(sampler, ddof=1) / np.sqrt(len(sampler))
    elif paired:
        diffs = x - y
        for i in range(R):
            sampler = np.random.choice(diffs, len(diffs), replace=True)
            m_vec[i] = np.mean(sampler)
            m_se_vec[i] = np.std(sampler, ddof=1) / np.sqrt(len(sampler))
    else:
        for i in range(R):
            x_boot = np.random.choice(x, len(x), replace=True)
            y_boot = np.random.choice(y, len(y), replace=True)
            m_vec[i] = np.mean(x_boot) - np.mean(y_boot)
            m_se_vec[i] = np.sqrt(np.var(x_boot, ddof=1)/len(x_boot) + np.var(y_boot, ddof=1)/len(y_boot))

    tstat = null_test['statistic']
    if y is None:
        TSTAT = (m_vec - np.mean(x)) / m_se_vec
    elif paired:
        TSTAT = (m_vec - np.mean(x-y)) / m_se_vec
    else:
        TSTAT = (m_vec - (np.mean(x) - np.mean(y))) / m_se_vec

    TSTAT = TSTAT[np.isfinite(TSTAT)]

    if alternative == 'two.sided':
        p_value = 2 * min(np.mean(TSTAT <= tstat), np.mean(TSTAT > tstat))
    elif alternative == 'greater':
        p_value = np.mean(TSTAT >= tstat)
    elif alternative == 'less':
        p_value = np.mean(TSTAT <= tstat)
    else:
        p_value = np.nan
        warnings.warn("Equivalence and MET alternatives not fully implemented for boot_t_test.")

    if boot_ci == "stud":
        ci = _stud_ci(m_vec, m_se_vec, null_test['stderr'], np.mean(x) if y is None else np.mean(x-y), alpha)
    elif boot_ci == "basic":
        ci = _basic_ci(m_vec, np.mean(x) if y is None else np.mean(x-y), alpha)
    else: # perc
        ci = _perc_ci(m_vec, alpha)

    return {
        'p.value': p_value,
        'conf.int': ci,
        'estimate': np.mean(x) if y is None else np.mean(x-y) if paired else np.mean(x) - np.mean(y),
        'null.value': mu,
        'alternative': alternative,
        'stderr': np.std(m_vec)
    }

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
