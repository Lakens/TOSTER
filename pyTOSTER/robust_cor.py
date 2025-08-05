import numpy as np
from scipy import stats

def winval(x, tr=0.2):
    n = len(x)
    ibot = int(tr * n)
    itop = n - ibot
    y = np.sort(x)
    xbot = y[ibot]
    xtop = y[itop-1]
    winval = np.where(x <= xbot, xbot, x)
    winval = np.where(winval >= xtop, xtop, winval)
    return winval

def wincor(x, y, tr=0.2):
    x = np.array(x)
    y = np.array(y)
    mask = ~np.isnan(x) & ~np.isnan(y)
    x = x[mask]
    y = y[mask]

    xvec = winval(x, tr)
    yvec = winval(y, tr)
    r, p = stats.pearsonr(xvec, yvec)
    return r

def pbos(x, beta=0.2):
    n = len(x)
    temp = np.sort(np.abs(x - np.median(x)))
    omhatx = temp[int((1 - beta) * n)]
    if omhatx == 0: return np.median(x)
    psi = (x - np.median(x)) / omhatx
    i1 = len(psi[psi < -1])
    i2 = len(psi[psi > 1])
    sx = np.where(psi < -1, 0, x)
    sx = np.where(psi > 1, 0, sx)
    pbos_val = (np.sum(sx) + omhatx * (i2 - i1)) / (n - i1 - i2)
    return pbos_val

def pbcor(x, y, beta=0.2):
    x = np.array(x)
    y = np.array(y)
    mask = ~np.isnan(x) & ~np.isnan(y)
    x = x[mask]
    y = y[mask]

    temp = np.sort(np.abs(x - np.median(x)))
    omhatx = temp[int((1 - beta) * len(x))]
    temp = np.sort(np.abs(y - np.median(y)))
    omhaty = temp[int((1 - beta) * len(y))]

    if omhatx == 0 or omhaty == 0: return 0.0

    a = (x - pbos(x, beta)) / omhatx
    b = (y - pbos(y, beta)) / omhaty
    a = np.clip(a, -1, 1)
    b = np.clip(b, -1, 1)

    num = np.sum(a * b)
    den = np.sqrt(np.sum(a**2) * np.sum(b**2))
    if den == 0: return 0.0

    return num / den
