import numpy as np
from scipy import stats
from .robust_cor import wincor, pbcor

def boot_cor_test(x, y,
                  alternative='two.sided',
                  method='pearson',
                  alpha=0.05,
                  null=0,
                  R=1999,
                  tr=0.2,
                  beta=0.2):
    """
    Bootstrapped Correlation Coefficients.
    """
    if not isinstance(x, (list, np.ndarray)) or not isinstance(y, (list, np.ndarray)):
        raise TypeError("x and y must be list or numpy array.")

    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")

    # Remove missing values
    x, y = np.array(x), np.array(y)

    mask = ~np.isnan(x) & ~np.isnan(y)
    x, y = x[mask], y[mask]

    n = len(x)

    # Bootstrap resampling
    bvec = np.zeros(R)

    for i in range(R):
        idx = np.random.choice(n, n, replace=True)

        if method == 'pearson':
            bvec[i] = stats.pearsonr(x[idx], y[idx])[0]
        elif method == 'spearman':
            bvec[i] = stats.spearmanr(x[idx], y[idx])[0]
        elif method == 'kendall':
            bvec[i] = stats.kendalltau(x[idx], y[idx])[0]
        elif method == 'winsorized':
            bvec[i] = wincor(x[idx], y[idx], tr=tr)
        elif method == 'bendpercent':
            bvec[i] = pbcor(x[idx], y[idx], beta=beta)
        else:
            raise ValueError(f"Method '{method}' not supported.")

    # Calculate p-value
    if alternative == 'two.sided':
        phat = (np.sum(bvec < null) + 0.5 * np.sum(bvec == null)) / R
        p_value = 2 * min(phat, 1 - phat)
    elif alternative == 'greater':
        p_value = np.sum(bvec <= null) / R
    elif alternative == 'less':
        p_value = np.sum(bvec >= null) / R
    elif alternative == 'equivalence':
        if isinstance(null, (int, float)):
            low_eqbound, high_eqbound = -abs(null), abs(null)
        else:
            low_eqbound, high_eqbound = min(null), max(null)
        p1 = np.sum(bvec >= high_eqbound) / R
        p2 = np.sum(bvec <= low_eqbound) / R
        p_value = max(p1, p2)
    elif alternative == 'minimal.effect':
        if isinstance(null, (int, float)):
            low_eqbound, high_eqbound = -abs(null), abs(null)
        else:
            low_eqbound, high_eqbound = min(null), max(null)
        p1 = np.sum(bvec <= high_eqbound) / R
        p2 = np.sum(bvec >= low_eqbound) / R
        p_value = min(p1, p2)
    else:
        raise ValueError(f"Alternative '{alternative}' not supported.")

    # Confidence interval
    if alternative in ['equivalence', 'minimal.effect']:
        ci_level = 1 - 2 * alpha
    else:
        ci_level = 1 - alpha

    conf_int = np.quantile(bvec, [(1 - ci_level) / 2, 1 - (1 - ci_level) / 2])

    # Initial correlation
    if method == 'pearson':
        estimate, _ = stats.pearsonr(x, y)
    elif method == 'spearman':
        estimate, _ = stats.spearmanr(x, y)
    elif method == 'kendall':
        estimate, _ = stats.kendalltau(x, y)
    elif method == 'winsorized':
        estimate = wincor(x, y, tr=tr)
    elif method == 'bendpercent':
        estimate = pbcor(x, y, beta=beta)

    stderr = np.std(bvec)

    results = {
        'p.value': p_value,
        'parameter': {'n': n},
        'conf.int': conf_int,
        'estimate': estimate,
        'stderr': stderr,
        'null.value': null,
        'alternative': alternative,
        'method': f"Bootstrapped {method} correlation",
        'boot_res': bvec
    }
    return results
