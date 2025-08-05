import numpy as np
import pandas as pd
from .gg_curv_t import gg_curv_t
from .cor_test import _cor_to_ci

def _corr_curv(r, n, method='pearson', steps=5000):
    intrvls = np.linspace(0, 1, steps + 1)
    intrvls = intrvls[(intrvls > 0) & (intrvls < 1)]

    results = [_cor_to_ci(r, n, ci=i, method=method) for i in intrvls]

    df = pd.DataFrame(results, columns=['lower.limit', 'upper.limit'])
    df['intrvl.level'] = intrvls
    df['pvalue'] = 1 - intrvls

    densdf = pd.DataFrame({'x': np.concatenate([df['lower.limit'], df['upper.limit']])})

    return [df, densdf]

def plot_cor(r, n, method='pearson', type=["c", "cd"], levels=[0.68, 0.9, 0.95, 0.999]):
    dat = _corr_curv(r, n, method)

    return gg_curv_t(dat, type=type, levels=levels,
                     xaxis="Correlation Coefficient",
                     yaxis1="p-value",
                     yaxis2="Confidence Interval (%)")
