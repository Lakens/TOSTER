import numpy as np
import pandas as pd
from .gg_curv_t import gg_curv_t
from .pes_calc import pes_ci

def _pes_curv(Fstat, df1, df2, steps=5000):
    intrvls = np.linspace(0, 1, steps + 1)
    intrvls = intrvls[(intrvls > 0) & (intrvls < 1)]

    results = [pes_ci(Fstat, df1, df2, conf_level=i) for i in intrvls]

    df = pd.DataFrame(results, columns=['lower.limit', 'upper.limit'])
    df['intrvl.level'] = intrvls
    df['pvalue'] = 1 - intrvls

    densdf = pd.DataFrame({'x': np.concatenate([df['lower.limit'], df['upper.limit']])})

    return [df, densdf]

def plot_pes(Fstat, df1, df2, type=["c", "cd"], levels=[0.68, 0.9, 0.95, 0.999]):
    dat = _pes_curv(Fstat, df1, df2)

    return gg_curv_t(dat, type=type, levels=levels,
                     xaxis="partial eta-squared",
                     yaxis1="p-value",
                     yaxis2="Confidence Interval (%)")
