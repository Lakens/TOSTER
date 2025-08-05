import numpy as np
from scipy import stats
import pandas as pd
from .cohend_calcs import d_est_one, d_est_pair, d_est_ind

def t_tost(x, y=None, hypothesis="EQU", paired=False, var_equal=False, eqb=None,
           low_eqbound=None, high_eqbound=None, eqbound_type="raw", alpha=0.05,
           bias_correction=True, rm_correction=False, glass=None, mu=0):

    if eqb is not None:
        if isinstance(eqb, (int, float)):
            high_eqbound = abs(eqb)
            low_eqbound = -abs(eqb)
        else:
            high_eqbound = max(eqb)
            low_eqbound = min(eqb)

    smd_type = 'g' if bias_correction else 'd'

    x = np.array(x)[~np.isnan(np.array(x))]
    if y is not None:
        y = np.array(y)[~np.isnan(np.array(y))]

    if y is None: # one-sample
        res = stats.ttest_1samp(x, popmean=mu)
        smd_res = d_est_one(len(x), np.mean(x) - mu, np.std(x, ddof=1), smd_type)
        if eqbound_type == 'SMD':
            low_eqbound_raw = low_eqbound * np.std(x, ddof=1)
            high_eqbound_raw = high_eqbound * np.std(x, ddof=1)
        else:
            low_eqbound_raw = low_eqbound
            high_eqbound_raw = high_eqbound

        t_low = stats.ttest_1samp(x, popmean=low_eqbound_raw, alternative='greater')
        t_high = stats.ttest_1samp(x, popmean=high_eqbound_raw, alternative='less')
        df = len(x) - 1

    elif paired:
        res = stats.ttest_rel(x, y)
        diffs = x - y
        r12, _ = stats.pearsonr(x, y)
        smd_res = d_est_pair(len(x), np.mean(x), np.mean(y), np.std(x, ddof=1), np.std(y, ddof=1), r12, smd_type, rm_correction)
        if eqbound_type == 'SMD':
            sd_diff = np.sqrt(np.std(x, ddof=1)**2 + np.std(y, ddof=1)**2 - 2 * r12 * np.std(x, ddof=1) * np.std(y, ddof=1))
            low_eqbound_raw = low_eqbound * sd_diff
            high_eqbound_raw = high_eqbound * sd_diff
        else:
            low_eqbound_raw = low_eqbound
            high_eqbound_raw = high_eqbound

        t_low = stats.ttest_1samp(diffs, popmean=low_eqbound_raw, alternative='greater')
        t_high = stats.ttest_1samp(diffs, popmean=high_eqbound_raw, alternative='less')
        df = len(x) - 1
    else: # independent
        res = stats.ttest_ind(x, y, equal_var=var_equal)
        smd_res = d_est_ind(len(x), len(y), np.mean(x), np.mean(y), np.std(x, ddof=1), np.std(y, ddof=1), smd_type, var_equal, glass)
        if eqbound_type == 'SMD':
            if var_equal:
                sd_pooled = np.sqrt((((len(x) - 1) * np.std(x, ddof=1)**2) + ((len(y) - 1) * np.std(y, ddof=1)**2)) / (len(x) + len(y) - 2))
            else:
                sd_pooled = np.sqrt((np.std(x, ddof=1)**2 + np.std(y, ddof=1)**2) / 2)
            low_eqbound_raw = low_eqbound * sd_pooled
            high_eqbound_raw = high_eqbound * sd_pooled
        else:
            low_eqbound_raw = low_eqbound
            high_eqbound_raw = high_eqbound

        t_low = stats.ttest_ind(x, y, equal_var=var_equal, alternative='greater')
        t_high = stats.ttest_ind(x, y, equal_var=var_equal, alternative='less')
        if var_equal:
            df = len(x) + len(y) - 2
        else:
            df = res.df

    if hypothesis == 'EQU':
        pTOST = max(t_low.pvalue, t_high.pvalue)
    else:
        pTOST = min(t_low.pvalue, t_high.pvalue)

    return {
        'TOST': pd.DataFrame({
            't': [res.statistic, t_low.statistic, t_high.statistic],
            'p.value': [res.pvalue, t_low.pvalue, t_high.pvalue],
            'df': [df, df, df]
        }, index=['t-test', 'TOST Lower', 'TOST Upper']),
        'smd': smd_res,
        'effsize': pd.DataFrame({
            'estimate': [res.statistic * np.sqrt(1/len(x) + 1/len(y)) if y is not None and not paired else np.mean(x-y) if paired else np.mean(x), smd_res['d']],
            'SE': [np.sqrt(1/len(x) + 1/len(y)) if y is not None and not paired else np.std(x-y, ddof=1)/np.sqrt(len(x)) if paired else np.std(x, ddof=1)/np.sqrt(len(x)), smd_res['se']]
        }, index=['Raw', smd_type]),
        'decision': {'TOST_p': pTOST}
    }
