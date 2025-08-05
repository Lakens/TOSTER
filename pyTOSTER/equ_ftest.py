import numpy as np
from scipy import stats
from scipy.optimize import brentq
import warnings

def _get_ncp_f(F, df1, df2, conf_level=0.95):
    alpha = 1 - conf_level

    # Lower limit
    func_lower = lambda ncp: stats.f.cdf(F, df1, df2, nc=ncp) - (1-alpha/2)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lower_ncp = brentq(func_lower, a=0, b=1000)
    except ValueError:
        lower_ncp = 0

    # Upper limit
    func_upper = lambda ncp: stats.f.cdf(F, df1, df2, nc=ncp) - (alpha/2)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            upper_ncp = brentq(func_upper, a=0, b=1000)
    except ValueError:
        upper_ncp = 1000

    return lower_ncp, upper_ncp

def equ_ftest(Fstat, df1, df2, eqbound, MET=False, alpha=0.05):

    pes = Fstat * df1 / (Fstat * df1 + df2)

    lower_ncp, upper_ncp = _get_ncp_f(Fstat, df1, df2, 1 - alpha)

    ll_pes = lower_ncp / (lower_ncp + df1 + df2 + 1)
    ul_pes = upper_ncp / (upper_ncp + df1 + df2 + 1)

    if np.isnan(ll_pes): ll_pes = 0
    if np.isnan(ul_pes): ul_pes = 1

    f2 = eqbound / (1 - eqbound)
    lambda_ = f2 * (df1 + df2 + 1)

    if MET:
        p_value = 1 - stats.f.cdf(Fstat, df1, df2, nc=lambda_)
        method = "Minimal Effect Test from F-test"
    else:
        p_value = stats.f.cdf(Fstat, df1, df2, nc=lambda_)
        method = "Equivalence Test from F-test"

    return {
        'statistic': Fstat,
        'p.value': p_value,
        'conf.int': [ll_pes, ul_pes],
        'estimate': pes,
        'null.value': eqbound,
        'alternative': None,
        'method': method
    }
