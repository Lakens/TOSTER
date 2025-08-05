import numpy as np
from scipy import stats

def extract_r_paired(m1, sd1, m2, sd2=None, n=None, tstat=None, pvalue=None):
    if n is None:
        raise ValueError("Sample size `n` must be provided.")
    if sd2 is None:
        sd2 = sd1
    if tstat is None:
        if pvalue is None:
            raise ValueError("tstat or pvalue must be provided")
        tstat = stats.t.ppf(1 - pvalue / 2, n - 1)

    corr = (sd2**2 + sd1**2 - tstat**(-2) * n * (m1 - m2)**2) / (2 * sd2 * sd1)
    return corr
