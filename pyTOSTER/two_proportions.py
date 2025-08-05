import numpy as np
from scipy import stats

def _p_from_z(z, alternative='two.sided'):
    if alternative == 'two.sided':
        return 2 * (1 - stats.norm.cdf(abs(z)))
    elif alternative == 'greater':
        return 1 - stats.norm.cdf(z)
    elif alternative == 'less':
        return stats.norm.cdf(z)
    else:
        raise ValueError("alternative must be 'two.sided', 'greater', or 'less'")

def test_prop_dif(p1, p2, n1, n2, null, alternative, alpha):
    if null is None: null = 0
    prop_dif = p1 - p2
    prop_se = np.sqrt((p1 * (1 - p1)) / n1 + (p2 * (1 - p2)) / n2)

    if alternative in ['equivalence', 'minimal.effect']:
        if isinstance(null, (int, float)):
            null = [-abs(null), abs(null)]

        lo_ztest = (prop_dif - min(null))/prop_se
        hi_ztest = (prop_dif - max(null))/prop_se

        if alternative == 'equivalence':
            lo_pvalue = _p_from_z(lo_ztest, alternative = "greater")
            hi_pvalue = _p_from_z(hi_ztest, alternative = "less")
            pval = max(lo_pvalue, hi_pvalue)
            ztest = hi_ztest if hi_pvalue >= lo_pvalue else lo_ztest
        else: # minimal.effect
            lo_pvalue = _p_from_z(lo_ztest, alternative = "less")
            hi_pvalue = _p_from_z(hi_ztest, alternative = "greater")
            pval = min(lo_pvalue, hi_pvalue)
            ztest = hi_ztest if hi_pvalue <= lo_pvalue else lo_ztest

        conf_level = 1 - alpha * 2
    else:
        ztest = (prop_dif - null) / prop_se
        pval = _p_from_z(ztest, alternative)
        conf_level = 1 - alpha if alternative == 'two.sided' else 1 - 2 * alpha

    z_mult = stats.norm.ppf(1 - (1 - conf_level) / 2)
    cint = prop_dif + np.array([-1, 1]) * z_mult * prop_se

    return {'STATISTIC': ztest, 'PVAL': pval, 'NVAL': null, 'ESTIMATE': prop_dif, 'CINT': cint, 'METHOD': 'z-test for difference in proportions'}

def twoprop_test(p1, p2, n1, n2, null=None, alpha=0.05, alternative='two.sided', effect_size='difference'):
    if effect_size != 'difference':
        raise NotImplementedError("Only 'difference' effect size is supported for twoprop_test.")

    res = test_prop_dif(p1, p2, n1, n2, null, alternative, alpha)

    return {
        'statistic': res['STATISTIC'],
        'p.value': res['PVAL'],
        'estimate': res['ESTIMATE'],
        'null.value': res['NVAL'],
        'conf.int': res['CINT'],
        'alternative': alternative,
        'method': res['METHOD']
    }
