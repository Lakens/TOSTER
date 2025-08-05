import numpy as np
from scipy import stats
import pandas as pd
import warnings
from .cohend_calcs import d_est_one, d_est_pair, d_est_ind

def _t_tost_helper(x, y=None, hypothesis="EQU", paired=False, var_equal=False, eqb=None, mu=0, bias_correction=True, rm_correction=False, glass=None):
    smd_type = 'g' if bias_correction else 'd'

    if len(eqb) == 1:
        low_eqbound, high_eqbound = -abs(eqb), abs(eqb)
    else:
        low_eqbound, high_eqbound = min(eqb), max(eqb)

    if y is None: # one-sample
        res = stats.ttest_1samp(x, popmean=mu)
        stderr = np.std(x, ddof=1) / np.sqrt(len(x))
        smd_res = d_est_one(len(x), np.mean(x) - mu, np.std(x, ddof=1), smd_type)
        t_low = stats.ttest_1samp(x, popmean=low_eqbound, alternative='greater').statistic
        t_high = stats.ttest_1samp(x, popmean=high_eqbound, alternative='less').statistic
        df = len(x) - 1
    elif paired:
        res = stats.ttest_rel(x, y)
        stderr = np.std(x - y, ddof=1) / np.sqrt(len(x))
        r12, _ = stats.pearsonr(x, y)
        smd_res = d_est_pair(len(x), np.mean(x), np.mean(y), np.std(x, ddof=1), np.std(y, ddof=1), r12, smd_type, rm_correction)
        t_low = stats.ttest_rel(x, y, alternative='greater').statistic
        t_high = stats.ttest_rel(x, y, alternative='less').statistic
        df = len(x) - 1
    else: # independent
        res = stats.ttest_ind(x, y, equal_var=var_equal)
        stderr = np.sqrt(np.var(x, ddof=1)/len(x) + np.var(y, ddof=1)/len(y))
        smd_res = d_est_ind(len(x), len(y), np.mean(x), np.mean(y), np.std(x, ddof=1), np.std(y, ddof=1), smd_type, var_equal, glass)
        t_low = stats.ttest_ind(x, y, equal_var=var_equal, alternative='greater').statistic
        t_high = stats.ttest_ind(x, y, equal_var=var_equal, alternative='less').statistic
        if var_equal:
            df = len(x) + len(y) - 2
        else:
            df = (np.var(x, ddof=1)/len(x) + np.var(y, ddof=1)/len(y))**2 / \
                 ((np.var(x, ddof=1)/len(x))**2/(len(x)-1) + (np.var(y, ddof=1)/len(y))**2/(len(y)-1))

    return {
        'smd': smd_res,
        'effsize': {'estimate': [res.statistic * stderr, smd_res['d']], 'SE': [stderr, smd_res['se']]},
        'TOST': {'t': [res.statistic, t_low, t_high], 'df': df}
    }


def boot_t_tost(x, y=None, hypothesis="EQU", paired=False, var_equal=False, eqb=None, alpha=0.05, mu=0, bias_correction=True, rm_correction=False, glass=None, boot_ci="stud", R=1999):

    nullTOST = _t_tost_helper(x, y, hypothesis, paired, var_equal, eqb, mu, bias_correction, rm_correction, glass)

    d_vec = np.zeros(R)
    m_vec = np.zeros(R)
    d_se_vec = np.zeros(R)
    m_se_vec = np.zeros(R)

    for i in range(R):
        if y is None:
            sampler = np.random.choice(x, len(x), replace=True)
            res = _t_tost_helper(sampler, None, hypothesis, paired, var_equal, eqb, mu, bias_correction, rm_correction, glass)
        elif paired:
            sampler = np.random.choice(len(x), len(x), replace=True)
            res = _t_tost_helper(x[sampler], y[sampler], hypothesis, paired, var_equal, eqb, mu, bias_correction, rm_correction, glass)
        else:
            x_boot = np.random.choice(x, len(x), replace=True)
            y_boot = np.random.choice(y, len(y), replace=True)
            res = _t_tost_helper(x_boot, y_boot, hypothesis, paired, var_equal, eqb, mu, bias_correction, rm_correction, glass)

        d_vec[i] = res['smd']['d']
        m_vec[i] = res['effsize']['estimate'][0]
        d_se_vec[i] = res['effsize']['SE'][1]
        m_se_vec[i] = res['effsize']['SE'][0]

    tstat, tstat_l, tstat_u = nullTOST['TOST']['t']

    if paired:
        TSTAT = (m_vec - np.mean(x-y)) / m_se_vec
    elif y is None:
        TSTAT = (m_vec - np.mean(x)) / m_se_vec
    else:
        TSTAT = (m_vec - (np.mean(x) - np.mean(y))) / m_se_vec

    TSTAT = TSTAT[np.isfinite(TSTAT)]

    boot_pval = 2 * min(np.mean(TSTAT <= tstat), np.mean(TSTAT > tstat))
    p_l = np.mean(TSTAT > tstat_l)
    p_u = np.mean(TSTAT < tstat_u)

    if hypothesis == "EQU":
        pTOST = max(p_l, p_u)
    else:
        pTOST = min(p_l, p_u)

    if boot_ci == "stud":
        boot_cint = _stud_ci(m_vec, m_se_vec, nullTOST['effsize']['SE'][0], nullTOST['effsize']['estimate'][0], alpha*2)
        d_cint = _stud_ci(d_vec, d_se_vec, nullTOST['effsize']['SE'][1], nullTOST['effsize']['estimate'][1], alpha*2)
    elif boot_ci == "basic":
        boot_cint = _basic_ci(m_vec, nullTOST['effsize']['estimate'][0], alpha*2)
        d_cint = _basic_ci(d_vec, nullTOST['effsize']['estimate'][1], alpha*2)
    else: # perc
        boot_cint = _perc_ci(m_vec, alpha*2)
        d_cint = _perc_ci(d_vec, alpha*2)

    return {
        'p.value': boot_pval,
        'TOST_p.value': pTOST,
        'conf.int': boot_cint,
        'smd_conf.int': d_cint,
        'estimate': nullTOST['effsize']['estimate'][0],
        'smd_estimate': nullTOST['effsize']['estimate'][1]
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
