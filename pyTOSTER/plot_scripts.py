import numpy as np
import pandas as pd
from scipy import stats

def _t_ci(tost_res, alpha):
    est = tost_res['effsize']['estimate'][0]
    se = tost_res['effsize']['SE'][0]
    df = tost_res['TOST']['df'][0]
    t_val = stats.t.ppf(1 - alpha / 2, df)
    return est - t_val * se, est + t_val * se

def _d_ci(d, df, ncp, sigma, smd_ci, t_stat, alpha):
    # Simplified version of d_CI
    if sigma is np.nan:
        return np.nan, np.nan
    return d - 1.96 * sigma, d + 1.96 * sigma

def t_curv(tost_res, steps=5000):
    intrvls = np.linspace(0, 1, steps + 1)
    intrvls = intrvls[(intrvls > 0) & (intrvls < 1)]

    results = [_t_ci(tost_res, 1-i) for i in intrvls]

    df = pd.DataFrame(results, columns=['lower.limit', 'upper.limit'])
    df['intrvl.level'] = intrvls
    df['pvalue'] = 1 - intrvls

    densdf = pd.DataFrame({'x': np.concatenate([df['lower.limit'], df['upper.limit']])})

    return [df, densdf]

def d_curv(tost_res, steps=5000):
    intrvls = np.linspace(0, 1, steps + 1)
    intrvls = intrvls[(intrvls > 0) & (intrvls < 1)]

    smd = tost_res['smd']
    results = [_d_ci(smd['d'], smd.get('d_df'), smd.get('d_lambda'), smd.get('d_sigma'), smd.get('smd_ci'), smd.get('t_stat'), 1-i) for i in intrvls]

    df = pd.DataFrame(results, columns=['lower.limit', 'upper.limit'])
    df['intrvl.level'] = intrvls
    df['pvalue'] = 1 - intrvls

    densdf = pd.DataFrame({'x': np.concatenate([df['lower.limit'], df['upper.limit']])})

    return [df, densdf]

def d_curv_raw(d, df, ncp, sigma, smd_ci="goulet", steps=5000):
    intrvls = np.linspace(0, 1, steps + 1)
    intrvls = intrvls[(intrvls > 0) & (intrvls < 1)]

    results = [_d_ci(d, df, ncp, sigma, smd_ci, None, 1-i) for i in intrvls]

    df = pd.DataFrame(results, columns=['lower.limit', 'upper.limit'])
    df['intrvl.level'] = intrvls
    df['pvalue'] = 1 - intrvls

    densdf = pd.DataFrame({'x': np.concatenate([df['lower.limit'], df['upper.limit']])})

    return [df, densdf]
