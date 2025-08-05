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

def compare_cor(r1, df1, r2, df2, method='fisher', alternative='two.sided', null=0):

    if alternative in ['equivalence', 'minimal.effect']:
        if isinstance(null, (int, float)):
            null = [-abs(null), abs(null)]

    if method == 'fisher':
        z1 = _rho_to_z(r1)
        z2 = _rho_to_z(r2)
        diff = z1 - z2
        znull = _rho_to_z(np.array(null))
        z_se = np.sqrt(1/(df1 - 1) + 1/(df2 - 1))

        if alternative == 'equivalence':
            zlo = diff - min(znull)
            plo = _p_from_z(zlo / z_se, 'greater')
            zhi = diff - max(znull)
            phi = _p_from_z(zhi / z_se, 'less')
            p_value = max(plo, phi)
            z = zhi if phi >= plo else zlo
        elif alternative == 'minimal.effect':
            zlo = diff - min(znull)
            plo = _p_from_z(zlo / z_se, 'less')
            zhi = diff - max(znull)
            phi = _p_from_z(zhi / z_se, 'greater')
            p_value = min(plo, phi)
            z = zhi if phi <= plo else zlo
        else:
            z_diff = diff - znull
            z = z_diff / z_se
            p_value = _p_from_z(z, alternative)

    elif method == 'kraatz':
        se = np.sqrt((1 - r1**2)**2 / df1 + (1 - r2**2)**2 / df2)
        diff = r1 - r2

        if alternative == 'equivalence':
            zlo = diff - min(null)
            plo = _p_from_z(zlo / se, 'greater')
            zhi = diff - max(null)
            phi = _p_from_z(zhi / se, 'less')
            p_value = max(plo, phi)
            z = zhi if phi >= plo else zlo
        elif alternative == 'minimal.effect':
            zlo = diff - min(null)
            plo = _p_from_z(zlo / se, 'less')
            zhi = diff - max(null)
            phi = _p_from_z(zhi / se, 'greater')
            p_value = min(plo, phi)
            z = zhi if phi <= plo else zlo
        else:
            diff_null = diff - null
            z = diff_null / se
            p_value = _p_from_z(z, alternative)
    else:
        raise ValueError("method must be 'fisher' or 'kraatz'")

    return {
        'statistic': z,
        'p.value': p_value,
        'estimate': r1 - r2,
        'null.value': null,
        'alternative': alternative,
        'method': f'Difference between two independent correlations ({method})'
    }
