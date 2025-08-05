import numpy as np
import pandas as pd
from .cohend_calcs import d_est_one, d_est_pair, d_est_ind
from .tsum_funcs import tsum_test

def tsum_tost(m1, sd1, n1, m2=None, sd2=None, n2=None, r12=None,
              hypothesis="EQU", paired=False, var_equal=False, eqb=None,
              low_eqbound=None, high_eqbound=None, mu=0, eqbound_type="raw",
              alpha=0.05, bias_correction=True, rm_correction=False, glass=None):

    if eqb is not None:
        if isinstance(eqb, (int, float)):
            high_eqbound = abs(eqb)
            low_eqbound = -abs(eqb)
        else:
            high_eqbound = max(eqb)
            low_eqbound = min(eqb)

    smd_type = 'g' if bias_correction else 'd'

    if m2 is None: # one-sample
        smd_res = d_est_one(n1, m1 - mu, sd1, smd_type)
        d_denom = sd1
    elif paired:
        smd_res = d_est_pair(n1, m1, m2, sd1, sd2, r12, smd_type, rm_correction)
        d_denom = np.sqrt(sd1**2 + sd2**2 - 2 * r12 * sd1 * sd2)
    else: # independent
        smd_res = d_est_ind(n1, n2, m1, m2, sd1, sd2, smd_type, var_equal, glass)
        if var_equal:
            d_denom = np.sqrt((((n1 - 1) * sd1**2) + ((n2 - 1) * sd2**2)) / (n1 + n2 - 2))
        else:
            d_denom = np.sqrt((sd1**2 + sd2**2) / 2)


    if eqbound_type == 'SMD':
        low_eqbound_raw = low_eqbound * d_denom
        high_eqbound_raw = high_eqbound * d_denom
    else:
        low_eqbound_raw = low_eqbound
        high_eqbound_raw = high_eqbound

    tresult = tsum_test(m1=m1, sd1=sd1, n1=n1, m2=m2, sd2=sd2, n2=n2, r12=r12,
                        paired=paired, var_equal=var_equal, mu=mu,
                        conf_level=1-alpha*2, alternative="two.sided")

    t_low = tsum_test(m1=m1, sd1=sd1, n1=n1, m2=m2, sd2=sd2, n2=n2, r12=r12,
                      paired=paired, var_equal=var_equal, mu=low_eqbound_raw,
                      alternative='greater')

    t_high = tsum_test(m1=m1, sd1=sd1, n1=n1, m2=m2, sd2=sd2, n2=n2, r12=r12,
                       paired=paired, var_equal=var_equal, mu=high_eqbound_raw,
                       alternative='less')

    if hypothesis == 'EQU':
        pTOST = max(t_low['p.value'], t_high['p.value'])
    else:
        pTOST = min(t_low['p.value'], t_high['p.value'])

    return {
        'TOST': pd.DataFrame({
            't': [tresult['statistic'], t_low['statistic'], t_high['statistic']],
            'p.value': [tresult['p.value'], t_low['p.value'], t_high['p.value']],
            'df': [tresult['df'], t_low['df'], t_high['df']]
        }, index=['t-test', 'TOST Lower', 'TOST Upper']),
        'smd': smd_res,
        'effsize': pd.DataFrame({
            'estimate': [tresult['estimate'], smd_res['d']],
            'SE': [tresult['stderr'], smd_res['se']]
        }, index=['Raw', smd_type]),
        'decision': {'TOST_p': pTOST}
    }
