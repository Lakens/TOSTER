import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def tostr(n, r, low_eqbound_r, high_eqbound_r, alpha=0.05, plot=True, verbose=True):
    """
    TOST function for a correlations.
    """
    if low_eqbound_r >= high_eqbound_r:
        print("Warning: The lower bound is equal to or larger than the upper bound.")
    if n < 3:
        raise ValueError("The sample size should be larger than 2.")
    if not (0 < alpha < 1):
        raise ValueError("The alpha level should be a positive value between 0 and 1.")
    if not (-1 <= r <= 1):
        raise ValueError("The correlation should be a value between -1 and 1.")

    # Fisher's z-transformation
    z_r = np.arctanh(r)
    z_low = np.arctanh(low_eqbound_r)
    z_high = np.arctanh(high_eqbound_r)

    se = 1 / np.sqrt(n - 3)

    z1 = (z_r - z_low) / se
    z2 = (z_r - z_high) / se

    p1 = 1 - stats.norm.cdf(z1)
    p2 = stats.norm.cdf(z2)
    ptost = max(p1, p2)

    # NHST
    if n - 2 > 0:
        t_nhst = r * np.sqrt(n - 2) / np.sqrt(1 - r**2)
        pttest = 2 * (1 - stats.t.cdf(abs(t_nhst), df=n-2))
    else:
        pttest = np.nan

    # CIs
    z_ci90_low = z_r - stats.norm.ppf(1 - alpha) * se
    z_ci90_high = z_r + stats.norm.ppf(1 - alpha) * se
    ci90_low = np.tanh(z_ci90_low)
    ci90_high = np.tanh(z_ci90_high)

    z_ci95_low = z_r - stats.norm.ppf(1 - alpha / 2) * se
    z_ci95_high = z_r + stats.norm.ppf(1 - alpha / 2) * se
    ci95_low = np.tanh(z_ci95_low)
    ci95_high = np.tanh(z_ci95_high)

    test_outcome = "significant" if pttest < alpha else "non-significant"
    tost_outcome = "significant" if ptost < alpha else "non-significant"

    if plot:
        plt.figure()
        plt.ylim(0, 1)
        plot_min = min(ci95_low, low_eqbound_r) - max(ci95_high - ci95_low, high_eqbound_r - low_eqbound_r) / 10
        plot_max = max(ci95_high, high_eqbound_r) + max(ci95_high - ci95_low, high_eqbound_r - low_eqbound_r) / 10
        plt.xlim(plot_min, plot_max)

        plt.yticks([])
        plt.xlabel("Correlation")
        title = (f"Equivalence bounds {low_eqbound_r:.3f} and {high_eqbound_r:.3f}\\n"
                 f"r = {r:.3f}\\n"
                 f"TOST: {100*(1-alpha*2):.0f}% CI [{ci90_low:.3f}; {ci90_high:.3f}] {tost_outcome}\\n"
                 f"NHST: {100*(1-alpha):.0f}% CI [{ci95_low:.3f}; {ci95_high:.3f}] {test_outcome}")
        plt.title(title.replace('\\n', '\n'))

        plt.plot(r, 0.5, 's', markersize=8, color='black')
        plt.axvline(high_eqbound_r, linestyle='--', color='black')
        plt.axvline(low_eqbound_r, linestyle='--', color='black')
        plt.axvline(0, linestyle='--', color='grey')

        plt.hlines(0.5, ci90_low, ci90_high, linewidth=3, color='black')
        plt.hlines(0.5, ci95_low, ci95_high, linewidth=1, color='black')

        plt.show()

    if verbose:
        print("TOST results:")
        print(f"p-value lower bound: {p1:.3f}")
        print(f"p-value upper bound: {p2:.3f}")
        print("\\n")
        print("Equivalence bounds (r):")
        print(f"low eqbound: {low_eqbound_r:.4f}\\nhigh eqbound: {high_eqbound_r:.4f}")
        print("\\n")
        print("TOST confidence interval:")
        print(f"lower bound {100*(1-alpha*2):.0f}% CI: {ci90_low:.3f}\\nupper bound {100*(1-alpha*2):.0f}% CI:  {ci90_high:.3f}")
        print("\\n")
        print("NHST confidence interval:")
        print(f"lower bound {100*(1-alpha):.0f}% CI: {ci95_low:.3f}\\nupper bound {100*(1-alpha):.0f}% CI:  {ci95_high:.3f}")
        print("\\n")
        print("Equivalence Test Result:")
        print(f"The equivalence test was {tost_outcome}, p = {ptost:.3f}, given equivalence bounds of {low_eqbound_r:.3f} and {high_eqbound_r:.3f} and an alpha of {alpha}.")
        print("\\n")
        print("Null Hypothesis Test Result:")
        print(f"The null hypothesis test was {test_outcome}, p = {pttest:.3f}, given an alpha of {alpha}.")

        if pttest <= alpha and ptost <= alpha:
            combined_outcome = ("NHST: reject null significance hypothesis that the correlation is equal to 0 \\n"
                                "TOST: reject null equivalence hypothesis")
        elif pttest < alpha and ptost > alpha:
            combined_outcome = ("NHST: reject null significance hypothesis that the correlation is equal to 0 \\n"
                                "TOST: don't reject null equivalence hypothesis")
        elif pttest > alpha and ptost <= alpha:
            combined_outcome = ("NHST: don't reject null significance hypothesis that the correlation is equal to 0 \\n"
                                "TOST: reject null equivalence hypothesis")
        else: # pttest > alpha and ptost > alpha
            combined_outcome = ("NHST: don't reject null significance hypothesis that the correlation is equal to 0 \\n"
                                  "TOST: don't reject null equivalence hypothesis")

        print("\\n")
        print(combined_outcome)

    results = {
        'r': r,
        'TOST_p1': p1,
        'TOST_p2': p2,
        'alpha': alpha,
        'low_eqbound_r': low_eqbound_r,
        'high_eqbound_r': high_eqbound_r,
        'LL_CI_TOST': ci90_low,
        'UL_CI_TOST': ci90_high,
        'LL_CI_TTEST': ci95_low,
        'UL_CI_TTEST': ci95_high,
        'NHST_p': pttest
    }
    return results
