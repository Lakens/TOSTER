import numpy as np
from scipy import stats
import warnings

def tsum_test(m1, sd1=None, n1=None, m2=None, sd2=None, n2=None, r12=None,
              paired=False, alternative='two.sided', mu=0, var_equal=False,
              conf_level=0.95):

    alpha = 1 - conf_level

    if m2 is None: # one-sample
        df = n1 - 1
        se = sd1 / np.sqrt(n1)
        t_stat = (m1 - mu) / se
        estimate = m1
    elif paired:
        if n1 != n2:
            warnings.warn("Unequal number of pairs; using smallest n")
            n1 = min(n1, n2)
        df = n1 - 1
        sd_diff = np.sqrt(sd1**2 + sd2**2 - 2 * r12 * sd1 * sd2)
        se = sd_diff / np.sqrt(n1)
        t_stat = ((m1 - m2) - mu) / se
        estimate = m1 - m2
    else: # independent
        if var_equal:
            df = n1 + n2 - 2
            s_pooled = np.sqrt((((n1 - 1) * sd1**2) + ((n2 - 1) * sd2**2)) / df)
            se = s_pooled * np.sqrt(1/n1 + 1/n2)
        else:
            se = np.sqrt(sd1**2 / n1 + sd2**2 / n2)
            df_num = (sd1**2/n1 + sd2**2/n2)**2
            df_den = (sd1**4 / (n1**2 * (n1 - 1))) + (sd2**4 / (n2**2 * (n2 - 1)))
            df = df_num / df_den
        t_stat = ((m1 - m2) - mu) / se
        estimate = m1 - m2

    if alternative == 'two.sided':
        p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df))
        t_crit = stats.t.ppf(1 - alpha / 2, df)
        lower_ci = estimate - t_crit * se
        upper_ci = estimate + t_crit * se
    elif alternative == 'greater':
        p_value = 1 - stats.t.cdf(t_stat, df)
        t_crit = stats.t.ppf(1-alpha, df)
        lower_ci = estimate - t_crit * se
        upper_ci = np.inf
    else: # less
        p_value = stats.t.cdf(t_stat, df)
        t_crit = stats.t.ppf(1-alpha, df)
        lower_ci = -np.inf
        upper_ci = estimate + t_crit * se

    return {
        'statistic': t_stat,
        'df': df,
        'p.value': p_value,
        'estimate': estimate,
        'stderr': se,
        'conf.int': [lower_ci, upper_ci]
    }
