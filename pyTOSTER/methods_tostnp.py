import pandas as pd
import numpy as np
from .htest_helpers import as_htest, describe_htest

def _rounder_stat(number, digits=3):
    if pd.isna(number) or not isinstance(number, (int, float)):
        return np.nan
    cutoff = 1 * 10**(-digits)
    if abs(number) < cutoff and number != 0:
        return f'{number:.{digits}e}'
    else:
        return round(number, digits)

def _printable_pval(pval, digits=3):
    cutoff = 1 * 10**(-digits)
    if pval < cutoff:
        return f'p < {cutoff}'
    else:
        return f'p = {round(pval, digits)}'

def describe_tost(x, digits=3):
    # Simplified version of describe_TOST
    htest = as_htest(x)
    return describe_htest(htest, digits=digits)

def print_tostnp(x, digits=4):
    effsize = x['effsize']
    tost = x['TOST']

    tost['p.value'] = [f"< 0.001" if p < 0.001 else round(p, 3) for p in tost['p.value']]
    effsize['CI'] = [f"[{round(l, digits)}, {round(u, digits)}]" for l, u in zip(effsize['lower.ci'], effsize['upper.ci'])]

    print(x['method'])
    print(x['decision']['TOST'])
    print(x['decision']['test'])
    print(x['decision']['combined'])
    print("\nTOST Results")
    print(tost)
    print("\nEffect Sizes")
    print(effsize[['estimate', 'CI', 'conf.level']])

def describe_tostnp(x, digits=3):
    return describe_tost(x, digits)
