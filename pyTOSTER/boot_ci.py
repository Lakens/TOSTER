import numpy as np

def basic_ci(boots_est, t0, alpha):
    conf = 1 - alpha
    qq = np.quantile(boots_est, [(1 - conf) / 2, (1 + conf) / 2])
    return 2 * t0 - qq[::-1]

def perc_ci(boots_est, alpha):
    return np.quantile(boots_est, [alpha/2, 1-alpha/2])

def stud_ci(boots_est, boots_se, se0, t0, alpha):
    conf = 1 - alpha
    z = (boots_est - t0) / boots_se
    z = z[np.isfinite(z)]
    qq = np.quantile(z, [(1 - conf) / 2, (1 + conf) / 2])
    return t0 - se0 * qq[::-1]
