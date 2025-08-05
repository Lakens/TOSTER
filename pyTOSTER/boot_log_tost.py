import numpy as np
from scipy import stats
import warnings

def _log_tost_helper(x, y, paired, var_equal, null, eqb):
    if paired:
        if len(x) != len(y):
            raise ValueError("x and y must have the same length for a paired test.")
        log_diff = np.log(x) - np.log(y)
        res = stats.ttest_1samp(log_diff, popmean=np.log(null))
        se = np.std(log_diff, ddof=1) / np.sqrt(len(log_diff)) if len(log_diff) > 1 else 0
        d = np.exp(np.mean(log_diff))
        df = len(log_diff) - 1
    else:
        log_x, log_y = np.log(x), np.log(y)
        res = stats.ttest_ind(log_x, log_y, equal_var=var_equal)
        se = np.sqrt(np.var(log_x, ddof=1)/len(log_x) + np.var(log_y, ddof=1)/len(log_y)) if len(log_x) > 1 and len(log_y) > 1 else 0
        d = np.exp(np.mean(log_x) - np.mean(log_y))
        if var_equal:
            df = len(log_x) + len(log_y) - 2
        else:
            df = (np.var(log_x)/len(log_x) + np.var(log_y)/len(log_y))**2 / \
                 ((np.var(log_x)/len(log_x))**2/(len(log_x)-1) + (np.var(log_y)/len(log_y))**2/(len(log_y)-1))

    low_eq, high_eq = (1/eqb, eqb) if eqb > 1 else (eqb, 1/eqb)

    if paired:
        t_low = stats.ttest_1samp(log_diff, popmean=np.log(low_eq), alternative='greater').statistic
        t_high = stats.ttest_1samp(log_diff, popmean=np.log(high_eq), alternative='less').statistic
    else:
        t_low = stats.ttest_ind(log_x, log_y, equal_var=var_equal, alternative='greater').statistic
        t_high = stats.ttest_ind(log_x, log_y, equal_var=var_equal, alternative='less').statistic

    return {
        'smd': {'d': d},
        'effsize': {'estimate': [np.mean(res.statistic * se + np.log(null))], 'SE': [se]},
        'TOST': {'t': [res.statistic, t_low, t_high], 'df': df}
    }

def _perc_ci(boots_est, alpha):
    return np.quantile(boots_est, [alpha/2, 1-alpha/2])

def _basic_ci(boots_est, t0, alpha):
    conf = 1 - alpha
    qq = np.quantile(boots_est, [(1 - conf) / 2, (1 + conf) / 2])
    return 2 * t0 - qq[::-1]

def _stud_ci(boots_est, boots_se, se0, t0, alpha):
    conf = 1 - alpha
    z = (boots_est - t0) / boots_se
    z = z[np.isfinite(z)]
    qq = np.quantile(z, [(1 - conf) / 2, (1 + conf) / 2])
    return t0 - se0 * qq[::-1]

def boot_log_tost(x, y=None, hypothesis="EQU", paired=False, var_equal=False, eqb=1.25, alpha=0.05, null=1, boot_ci="stud", R=1999):
    if y is None:
        raise ValueError("One sample tests not supported.")
    if any(x <= 0) or any(y <= 0):
        raise ValueError("All values must be positive.")

    if len(eqb) == 1:
        high_eqbound, low_eqbound = (eqb, 1/eqb) if eqb > 1 else (1/eqb, eqb)
    else:
        high_eqbound, low_eqbound = max(eqb), min(eqb)

    log_x, log_y = np.log(x), np.log(y)

    if paired:
        if len(log_x) != len(log_y): raise ValueError("x and y must have the same length.")
        log_diff = log_x - log_y
        nullTOST = _log_tost_helper(x, y, paired, var_equal, null, eqb)
    else:
        nullTOST = _log_tost_helper(x, y, paired, var_equal, null, eqb)

    d_vec = np.zeros(R)
    m_vec = np.zeros(R)
    se_vec = np.zeros(R)

    for i in range(R):
        if paired:
            idx = np.random.choice(len(log_diff), len(log_diff), replace=True)
            res = _log_tost_helper(np.exp(log_x[idx]), np.exp(log_y[idx]), True, False, null, eqb)
        else:
            idx1 = np.random.choice(len(log_x), len(log_x), replace=True)
            idx2 = np.random.choice(len(log_y), len(log_y), replace=True)
            res = _log_tost_helper(np.exp(log_x[idx1]), np.exp(log_y[idx2]), False, var_equal, null, eqb)

        d_vec[i] = res['smd']['d']
        m_vec[i] = res['effsize']['estimate'][0]
        se_vec[i] = res['effsize']['SE'][0]

    # P-values
    tstat, tstat_l, tstat_u = nullTOST['TOST']['t']

    if paired:
        TSTAT = (m_vec - np.mean(log_diff)) / se_vec
    else:
        TSTAT = (m_vec - (np.mean(log_x) - np.mean(log_y))) / se_vec

    TSTAT = TSTAT[np.isfinite(TSTAT)]

    boot_pval = 2 * min(np.mean(TSTAT <= tstat), np.mean(TSTAT > tstat))
    p_l = np.mean(TSTAT > tstat_l)
    p_u = np.mean(TSTAT < tstat_u)

    if hypothesis == "EQU":
        pTOST = max(p_l, p_u)
    else:
        pTOST = min(p_l, p_u)

    # CIs
    if boot_ci == "stud":
        boot.cint = _stud_ci(m_vec, se_vec, nullTOST['effsize']['SE'][0], nullTOST['effsize']['estimate'][0], alpha*2)
    elif boot_ci == "basic":
        boot.cint = _basic_ci(m_vec, nullTOST['effsize']['estimate'][0], alpha*2)
    else: # perc
        boot.cint = _perc_ci(m_vec, alpha*2)

    d_cint = np.exp(boot.cint)

    return {
        'p.value': boot_pval,
        'TOST_p.value': pTOST,
        'conf.int': d_cint,
        'estimate': np.exp(nullTOST['effsize']['estimate'][0]),
        'boot_res': {'SMD': d_vec, 'raw': m_vec}
    }
