import pandas as pd
import numpy as np

def rounder_stat(number, digits=3):
    if pd.isna(number) or not isinstance(number, (int, float)):
        return np.nan
    cutoff = 1 * 10**(-digits)
    if abs(number) < cutoff and number != 0:
        return f'{number:.{digits}e}'
    else:
        return round(number, digits)

def printable_pval(pval, digits=3):
    cutoff = 1 * 10**(-digits)
    if pval < cutoff:
        return f'p < {cutoff}'
    else:
        return f'p = {round(pval, digits)}'

def df_htest(htest, test_statistics=True, show_ci=True):
    if not isinstance(htest, dict):
        raise TypeError("htest must be a dictionary.")

    df = pd.DataFrame({'method': [htest.get('method', '')]})

    if test_statistics:
        if htest.get('statistic') is not None:
            df['statistic'] = htest['statistic']
        if htest.get('parameter') is not None:
            df['parameter'] = htest['parameter']
        if htest.get('p.value') is not None:
            df['p.value'] = htest['p.value']

    if htest.get('estimate') is not None:
        df['estimate'] = htest['estimate']

    if show_ci and htest.get('conf.int') is not None:
        df['lower.ci'] = htest['conf.int'][0]
        df['upper.ci'] = htest['conf.int'][1]

    return df

def describe_htest(htest, alpha=None, digits=3):
    if not isinstance(htest, dict):
        raise TypeError("htest must be a dictionary.")

    if alpha is None:
        alpha = 0.05

    sig_state = "statistically significant" if htest['p.value'] < alpha else "not statistically significant"
    hyp_state = "The null hypothesis can be rejected." if htest['p.value'] < alpha else "The null hypothesis cannot be rejected."

    stat_state = ""
    if htest.get('statistic') is not None:
        stat_name = 'statistic'
        stat_est = htest['statistic']
        par_state = f"({htest['parameter']})" if htest.get('parameter') is not None else ""
        pval_state = printable_pval(htest['p.value'], digits)
        stat_state = f"{stat_name}{par_state} = {rounder_stat(stat_est, digits)}, {pval_state}"

    est_state = ""
    if htest.get('estimate') is not None and htest.get('conf.int') is not None:
        est_state = f"estimate = {rounder_stat(htest['estimate'], digits)}, CI [{rounder_stat(htest['conf.int'][0], digits)}, {rounder_stat(htest['conf.int'][1], digits)}]"

    return f"The {htest.get('method', '')} is {sig_state}, {stat_state}. {hyp_state} {est_state}"
