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

def test_odds_ratio(p1, p2, n1, n2, null, alternative, alpha):
    if null is None: null = 1
    q1 = 1 - p1
    q2 = 1 - p2
    OR = (p1 / q1) / (p2 / q2)
    se_logodds = np.sqrt(1/(n1*p1+0.5) + 1/(n1*q1+0.5) + 1/(n2*p2+0.5) + 1/(n2*q2+0.5))

    if alternative in ['equivalence', 'minimal.effect']:
        if isinstance(null, (int, float)):
            null = [null, 1/null]

        lo_ztest = (np.log(OR) - min(np.log(null))) / se_logodds
        hi_ztest = (np.log(OR) - max(np.log(null))) / se_logodds

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
        ztest = (np.log(OR) - np.log(null)) / se_logodds
        pval = _p_from_z(ztest, alternative)
        conf_level = 1 - alpha if alternative == 'two.sided' else 1 - 2 * alpha

    z_mult = stats.norm.ppf(1 - (1 - conf_level) / 2)
    cint = np.exp(np.log(OR) + np.array([-1, 1]) * z_mult * se_logodds)

    return {'STATISTIC': ztest, 'PVAL': pval, 'NVAL': null, 'ESTIMATE': OR, 'CINT': cint, 'METHOD': 'z-test for odds ratio'}

def test_risk_ratio(p1, p2, n1, n2, null, alternative, alpha):
    if null is None: null = 1
    phi = p1 / p2
    se_logrisk = np.sqrt((1-p1)/(n1*p1) + (1-p2)/(n2*p2))

    if alternative in ['equivalence', 'minimal.effect']:
        if isinstance(null, (int, float)):
            null = [null, 1/null]

        lo_ztest = (np.log(phi) - min(np.log(null))) / se_logrisk
        hi_ztest = (np.log(phi) - max(np.log(null))) / se_logrisk

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
        ztest = (np.log(phi) - np.log(null)) / se_logrisk
        pval = _p_from_z(ztest, alternative)
        conf_level = 1 - alpha if alternative == 'two.sided' else 1 - 2 * alpha

    z_mult = stats.norm.ppf(1 - (1 - conf_level) / 2)
    cint = phi * np.exp(np.array([-1, 1]) * z_mult * se_logrisk)

    return {'STATISTIC': ztest, 'PVAL': pval, 'NVAL': null, 'ESTIMATE': phi, 'CINT': cint, 'METHOD': 'z-test for risk ratio'}

def twoprop_test(p1, p2, n1, n2, null=None, alpha=0.05, alternative='two.sided', effect_size='difference'):
    if effect_size == 'difference':
        res = test_prop_dif(p1, p2, n1, n2, null, alternative, alpha)
    elif effect_size == 'odds.ratio':
        res = test_odds_ratio(p1, p2, n1, n2, null, alternative, alpha)
    elif effect_size == 'risk.ratio':
        res = test_risk_ratio(p1, p2, n1, n2, null, alternative, alpha)
    else:
        raise ValueError("effect_size must be 'difference', 'odds.ratio', or 'risk.ratio'")

    return {
        'statistic': res['STATISTIC'],
        'p.value': res['PVAL'],
        'estimate': res['ESTIMATE'],
        'null.value': res['NVAL'],
        'conf.int': res['CINT'],
        'alternative': alternative,
        'method': res['METHOD']
    }
