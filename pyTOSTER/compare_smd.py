import numpy as np
from scipy import stats
from .cohend_calcs import se_dz, se_ds

def _p_from_z(z, alternative='two.sided'):
    if alternative == 'two.sided':
        return 2 * (1 - stats.norm.cdf(abs(z)))
    elif alternative == 'greater':
        return 1 - stats.norm.cdf(z)
    elif alternative == 'less':
        return stats.norm.cdf(z)
    else:
        raise ValueError("alternative must be 'two.sided', 'greater', or 'less'")

def compare_smd(smd1, n1, smd2, n2, se1=None, se2=None, paired=False, alternative='two.sided', null=0):

    if alternative in ['equivalence', 'minimal.effect']:
        if isinstance(null, (int, float)):
            null = [-abs(null), abs(null)]

    if se1 is None:
        se1 = se_dz(smd1, n1) if paired else se_ds(smd1, n1)
    if se2 is None:
        se2 = se_dz(smd2, n2) if paired else se_ds(smd2, n2)

    se_diff = np.sqrt(se1**2 + se2**2)
    diff = smd1 - smd2

    if alternative == 'equivalence':
        null_low, null_high = min(null), max(null)
        z_low = (diff - null_low) / se_diff
        p_low = _p_from_z(z_low, 'greater')
        z_high = (diff - null_high) / se_diff
        p_high = _p_from_z(z_high, 'less')
        p_value = max(p_low, p_high)
        z = z_high if p_high >= p_low else z_low
    elif alternative == 'minimal.effect':
        null_low, null_high = min(null), max(null)
        z_low = (diff - null_low) / se_diff
        p_low = _p_from_z(z_low, 'less')
        z_high = (diff - null_high) / se_diff
        p_high = _p_from_z(z_high, 'greater')
        p_value = min(p_low, p_high)
        z = z_high if p_high <= p_low else z_low
    else:
        z = (diff - null) / se_diff
        p_value = _p_from_z(z, alternative)

    return {
        'statistic': z,
        'p.value': p_value,
        'estimate': diff,
        'null.value': null,
        'alternative': alternative,
        'method': 'Difference in SMDs'
    }
