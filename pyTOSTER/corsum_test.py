import numpy as np
from scipy import stats

def _rho_to_z(rho):
    return np.arctanh(rho)

def _p_from_z(z, alternative='two.sided'):
    if alternative == 'two.sided':
        return 2 * (1 - stats.norm.cdf(abs(z)))
    elif alternative == 'greater':
        return 1 - stats.norm.cdf(z)
    elif alternative == 'less':
        return stats.norm.cdf(z)
    else:
        raise ValueError("alternative must be 'two.sided', 'greater', or 'less'")

def _cor_to_ci(cor, n, ci=0.95, method='pearson'):
    z = _rho_to_z(cor)
    if method == 'pearson':
        se = 1 / np.sqrt(n - 3)
    elif method == 'spearman':
        se = 1.06 / np.sqrt(n - 3)
    elif method == 'kendall':
        se = 0.437 / np.sqrt(n - 4)

    z_ci = z + stats.norm.ppf([ (1-ci)/2, 1-(1-ci)/2 ]) * se
    return np.tanh(z_ci)

def corsum_test(r, n, alternative='two.sided', method='pearson', null=0, alpha=0.05):

    if alternative in ['equivalence', 'minimal.effect']:
        if isinstance(null, (int, float)):
            null = [-abs(null), abs(null)]

    z_r = _rho_to_z(r)
    z_null = _rho_to_z(np.array(null))

    if method == 'pearson':
        z_se = 1 / np.sqrt(n - 3)
    elif method == 'spearman':
        z_se = 1.06 / np.sqrt(n - 3)
    elif method == 'kendall':
        z_se = 0.437 / np.sqrt(n - 4)
    else:
        raise ValueError("method must be 'pearson', 'kendall', or 'spearman'")

    if alternative == 'equivalence':
        z_low = z_r - min(z_null)
        p_low = _p_from_z(z_low / z_se, 'greater')
        z_high = z_r - max(z_null)
        p_high = _p_from_z(z_high / z_se, 'less')
        p_value = max(p_low, p_high)
        z = z_high if p_high >= p_low else z_low
    elif alternative == 'minimal.effect':
        z_low = z_r - min(z_null)
        p_low = _p_from_z(z_low / z_se, 'less')
        z_high = z_r - max(z_null)
        p_high = _p_from_z(z_high / z_se, 'greater')
        p_value = min(p_low, p_high)
        z = z_high if p_high <= p_low else z_low
    else:
        z = z_r - z_null
        p_value = _p_from_z(z / z_se, alternative)

    ci = _cor_to_ci(r, n, 1-alpha*2 if alternative in ['equivalence', 'minimal.effect'] else 1-alpha, method)

    return {
        'statistic': z / z_se,
        'p.value': p_value,
        'conf.int': ci,
        'estimate': r,
        'null.value': null,
        'alternative': alternative,
        'method': f"Z-test for {method} correlation from summary"
    }
