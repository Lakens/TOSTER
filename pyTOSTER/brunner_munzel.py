import numpy as np
from scipy import stats

def brunner_munzel(x, y, alternative='two.sided', mu=0.5, paired=False, perm=False):
    """
    Brunner-Munzel test for stochastic equality of two samples.
    This is a simplified version of the R code. It does not support permutation tests.
    """
    if alternative not in ['two.sided', 'less', 'greater']:
        raise ValueError("alternative must be 'two.sided', 'less', or 'greater'")

    if perm:
        warnings.warn("Permutation test is not supported in this Python version.")

    x = np.array(x)[~np.isnan(np.array(x))]
    y = np.array(y)[~np.isnan(np.array(y))]

    n1 = len(x)
    n2 = len(y)

    if n1 == 0 or n2 == 0:
        return {'statistic': np.nan, 'p.value': np.nan, 'estimate': np.nan, 'df': np.nan}

    if paired:
        if n1 != n2:
            raise ValueError("x and y must have the same length for a paired test.")

        all_data = np.concatenate([y, x])
        N = len(all_data)

        rx = stats.rankdata(all_data)
        rix1 = stats.rankdata(y)
        rix2 = stats.rankdata(x)

        BM1 = (rx[:n1] - rix1) / n1
        BM2 = (rx[n1:] - rix2) / n2
        BM3 = BM1 - BM2

        pd = np.mean(BM2)
        m = np.mean(BM3)
        v = np.sum((BM3 - m)**2) / (n1 - 1)
        if v == 0: v = 1 / n1

        test_stat = np.sqrt(n1) * (pd - mu) / np.sqrt(v)
        df = n1 - 1
    else:
        rxy = stats.rankdata(np.concatenate([x, y]))
        rx = stats.rankdata(x)
        ry = stats.rankdata(y)

        pl2 = (rxy[:n1] - rx) / n2
        pl1 = (rxy[n1:] - ry) / n1

        pd = np.mean(pl2)

        s1 = np.var(pl2, ddof=1) / n1 if n1 > 1 else 0
        s2 = np.var(pl1, ddof=1) / n2 if n2 > 1 else 0

        V = (n1 + n2) * (s1 + s2)
        if V == 0: V = (n1 + n2) / (2 * n1 * n2)

        test_stat = np.sqrt(n1 + n2) * (pd - mu) / np.sqrt(V)
        df = (s1 + s2)**2 / (s1**2 / (n1 - 1) + s2**2 / (n2 - 1)) if n1 > 1 and n2 > 1 else 0

    if alternative == 'two.sided':
        p_value = 2 * min(stats.t.cdf(test_stat, df), 1 - stats.t.cdf(test_stat, df))
    elif alternative == 'greater':
        p_value = 1 - stats.t.cdf(test_stat, df)
    else: # less
        p_value = stats.t.cdf(test_stat, df)

    return {
        'statistic': test_stat,
        'p.value': p_value,
        'estimate': pd,
        'df': df
    }
