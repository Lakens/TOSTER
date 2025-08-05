import numpy as np
from scipy import stats
import warnings

try:
    import pingouin as pg
except ImportError:
    warnings.warn("pingouin is not installed. Please install it to use robust correlation methods.")
    pg = None

def boot_compare_cor(x1, y1, x2, y2,
                     alternative='two.sided',
                     method='pearson',
                     alpha=0.05,
                     null=0,
                     R=1999):
    """
    Comparing Correlations Between Independent Studies with Bootstrapping.

    Note: Robust methods 'winsorized' and 'bendpercent' require the 'pingouin' package.
    """
    if not isinstance(x1, (list, np.ndarray)) or not isinstance(y1, (list, np.ndarray)) or \
       not isinstance(x2, (list, np.ndarray)) or not isinstance(y2, (list, np.ndarray)):
        raise TypeError("x1, y1, x2, y2 must be list or numpy array.")

    if len(x1) != len(y1) or len(x2) != len(y2):
        raise ValueError("x and y must have the same length.")

    # Remove missing values
    x1, y1 = np.array(x1), np.array(y1)
    x2, y2 = np.array(x2), np.array(y2)

    mask1 = ~np.isnan(x1) & ~np.isnan(y1)
    x1, y1 = x1[mask1], y1[mask1]

    mask2 = ~np.isnan(x2) & ~np.isnan(y2)
    x2, y2 = x2[mask2], y2[mask2]

    n1, n2 = len(x1), len(x2)

    # Bootstrap resampling
    bvec1 = np.zeros(R)
    bvec2 = np.zeros(R)

    for i in range(R):
        idx1 = np.random.choice(n1, n1, replace=True)
        idx2 = np.random.choice(n2, n2, replace=True)

        if method == 'pearson':
            bvec1[i] = stats.pearsonr(x1[idx1], y1[idx1])[0]
            bvec2[i] = stats.pearsonr(x2[idx2], y2[idx2])[0]
        elif method == 'spearman':
            bvec1[i] = stats.spearmanr(x1[idx1], y1[idx1])[0]
            bvec2[i] = stats.spearmanr(x2[idx2], y2[idx2])[0]
        elif method == 'kendall':
            bvec1[i] = stats.kendalltau(x1[idx1], y1[idx1])[0]
            bvec2[i] = stats.kendalltau(x2[idx2], y2[idx2])[0]
        elif method == 'winsorized':
            if pg is None: raise ImportError("pingouin package is required for winsorized correlation.")
            bvec1[i] = pg.corr(stats.mstats.winsorize(x1[idx1]), stats.mstats.winsorize(y1[idx1]))['r'].values[0]
            bvec2[i] = pg.corr(stats.mstats.winsorize(x2[idx2]), stats.mstats.winsorize(y2[idx2]))['r'].values[0]
        elif method == 'bendpercent':
            if pg is None: raise ImportError("pingouin package is required for percentage bend correlation.")
            bvec1[i] = pg.corr(x1[idx1], y1[idx1], method='percbend')['r'].values[0]
            bvec2[i] = pg.corr(x2[idx2], y2[idx2], method='percbend')['r'].values[0]
        else:
            raise ValueError(f"Method '{method}' not supported.")

    bvec_diff = bvec1 - bvec2

    # Calculate p-value
    if alternative == 'two.sided':
        phat = (np.sum(bvec_diff < null) + 0.5 * np.sum(bvec_diff == null)) / R
        p_value = 2 * min(phat, 1 - phat)
    elif alternative == 'greater':
        p_value = np.sum(bvec_diff <= null) / R
    elif alternative == 'less':
        p_value = np.sum(bvec_diff >= null) / R
    elif alternative == 'equivalence':
        if isinstance(null, (int, float)):
            low_eqbound, high_eqbound = -abs(null), abs(null)
        else:
            low_eqbound, high_eqbound = min(null), max(null)
        p1 = np.sum(bvec_diff >= high_eqbound) / R
        p2 = np.sum(bvec_diff <= low_eqbound) / R
        p_value = max(p1, p2)
    elif alternative == 'minimal.effect':
        if isinstance(null, (int, float)):
            low_eqbound, high_eqbound = -abs(null), abs(null)
        else:
            low_eqbound, high_eqbound = min(null), max(null)
        p1 = np.sum(bvec_diff <= high_eqbound) / R
        p2 = np.sum(bvec_diff >= low_eqbound) / R
        p_value = min(p1, p2)
    else:
        raise ValueError(f"Alternative '{alternative}' not supported.")

    # Confidence interval
    if alternative in ['equivalence', 'minimal.effect']:
        ci_level = 1 - 2 * alpha
    else:
        ci_level = 1 - alpha

    conf_int = np.quantile(bvec_diff, [(1 - ci_level) / 2, 1 - (1 - ci_level) / 2])

    # Initial correlations
    if method == 'pearson':
        r1 = stats.pearsonr(x1, y1)[0]
        r2 = stats.pearsonr(x2, y2)[0]
    elif method == 'spearman':
        r1 = stats.spearmanr(x1, y1)[0]
        r2 = stats.spearmanr(x2, y2)[0]
    elif method == 'kendall':
        r1 = stats.kendalltau(x1, y1)[0]
        r2 = stats.kendalltau(x2, y2)[0]
    elif method == 'winsorized':
        if pg is None: raise ImportError("pingouin package is required for winsorized correlation.")
        r1 = pg.corr(stats.mstats.winsorize(x1), stats.mstats.winsorize(y1))['r'].values[0]
        r2 = pg.corr(stats.mstats.winsorize(x2), stats.mstats.winsorize(y2))['r'].values[0]
    elif method == 'bendpercent':
        if pg is None: raise ImportError("pingouin package is required for percentage bend correlation.")
        r1 = pg.corr(x1, y1, method='percbend')['r'].values[0]
        r2 = pg.corr(x2, y2, method='percbend')['r'].values[0]

    estimate = r1 - r2
    stderr = np.std(bvec_diff)

    results = {
        'p.value': p_value,
        'parameter': {'n1': n1, 'n2': n2},
        'conf.int': conf_int,
        'estimate': estimate,
        'stderr': stderr,
        'null.value': null,
        'alternative': alternative,
        'method': f"Bootstrapped difference in {method} correlation",
        'boot_res': {'diff': bvec_diff, 'r1': bvec1, 'r2': bvec2}
    }
    return results
