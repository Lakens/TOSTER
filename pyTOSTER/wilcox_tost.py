import numpy as np
from scipy import stats
import pandas as pd
from .rbs import np_ses

def wilcox_tost(x, y=None, hypothesis="EQU", paired=False, eqb=None,
                low_eqbound=None, high_eqbound=None, ses="rb", alpha=0.05, mu=0):

    if eqb is not None:
        if isinstance(eqb, (int, float)):
            high_eqbound = abs(eqb)
            low_eqbound = -abs(eqb)
        else:
            high_eqbound = max(eqb)
            low_eqbound = min(eqb)

    x = np.array(x)[~np.isnan(np.array(x))]
    if y is not None:
        y = np.array(y)[~np.isnan(np.array(y))]

    if y is None: # one-sample
        res = stats.wilcoxon(x - mu, alternative='two-sided')
        t_low = stats.wilcoxon(x - low_eqbound, alternative='greater')
        t_high = stats.wilcoxon(x - high_eqbound, alternative='less')
    elif paired:
        res = stats.wilcoxon(x, y, alternative='two-sided')
        t_low = stats.wilcoxon(x - y, x_alt=low_eqbound, alternative='greater')
        t_high = stats.wilcoxon(x - y, x_alt=high_eqbound, alternative='less')
    else: # independent
        res = stats.mannwhitneyu(x, y, alternative='two-sided')
        t_low = stats.mannwhitneyu(x - low_eqbound, y, alternative='greater')
        t_high = stats.mannwhitneyu(x - high_eqbound, y, alternative='less')

    if hypothesis == 'EQU':
        pTOST = max(t_low.pvalue, t_high.pvalue)
    else:
        pTOST = min(t_low.pvalue, t_high.pvalue)

    effsize = np_ses(x, y, mu, 1-alpha*2, paired, ses)

    return {
        'TOST': pd.DataFrame({
            'statistic': [res.statistic, t_low.statistic, t_high.statistic],
            'p.value': [res.pvalue, t_low.pvalue, t_high.pvalue]
        }, index=['NHST', 'TOST Lower', 'TOST Upper']),
        'effsize': effsize,
        'decision': {'TOST_p': pTOST}
    }
