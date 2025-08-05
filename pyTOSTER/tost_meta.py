import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def tost_meta(es, var=None, se=None, low_eqbound_d=None, high_eqbound_d=None, alpha=0.05, plot=True, verbose=True):
    """
    TOST function for meta-analysis.

    Parameters
    ----------
    es : float
        Meta-analytic effect size.
    var : float, optional
        Meta-analytic variance.
    se : float, optional
        Standard error.
    low_eqbound_d : float
        Lower equivalence bound in Cohen's d.
    high_eqbound_d : float
        Upper equivalence bound in Cohen's d.
    alpha : float, optional
        Alpha level (default is 0.05).
    plot : bool, optional
        Set whether results should be plotted (default is True).
    verbose : bool, optional
        Set whether text output should be generated (default is True).

    Returns
    -------
    dict
        A dictionary containing the TOST results.
    """
    if low_eqbound_d is None or high_eqbound_d is None:
        raise ValueError("Equivalence bounds must be specified.")

    if var is None and se is None:
        raise ValueError("Need to specify variance (var) or standard error (se).")

    if se is None:
        se = np.sqrt(var)

    if low_eqbound_d >= high_eqbound_d:
        print("Warning: The lower bound is equal to or larger than the upper bound.")

    if not (0 < alpha < 1):
        raise ValueError("The alpha level should be a positive value between 0 and 1.")

    z1 = (es - low_eqbound_d) / se
    p1 = 1 - stats.norm.cdf(z1)
    z2 = (es - high_eqbound_d) / se
    p2 = stats.norm.cdf(z2)

    z = es / se
    pttest = 2 * (1 - stats.norm.cdf(abs(z)))

    ll90 = es - stats.norm.ppf(1 - alpha) * se
    ul90 = es + stats.norm.ppf(1 - alpha) * se
    ll95 = es - stats.norm.ppf(1 - alpha / 2) * se
    ul95 = es + stats.norm.ppf(1 - alpha / 2) * se

    ptost = max(p1, p2)
    ztost = z1 if abs(z1) < abs(z2) else z2

    test_outcome = "significant" if pttest < alpha else "non-significant"
    tost_outcome = "significant" if ptost < alpha else "non-significant"

    if plot:
        plt.figure()
        plt.ylim(0, 1)
        plot_min = min(ll95, low_eqbound_d, es) - max(ul95 - ll95, high_eqbound_d - low_eqbound_d, es) / 10
        plot_max = max(ul95, high_eqbound_d, es) + max(ul95 - ll95, high_eqbound_d - low_eqbound_d, es) / 10
        plt.xlim(plot_min, plot_max)

        plt.yticks([])
        plt.xlabel("Effect size")
        title = (f"Equivalence bounds {low_eqbound_d:.3f} and {high_eqbound_d:.3f}\\n"
                 f"Effect size = {es:.3f}\\n"
                 f"TOST: {100*(1-alpha*2):.0f}% CI [{ll90:.3f}; {ul90:.3f}] {tost_outcome}\\n"
                 f"NHST: {100*(1-alpha):.0f}% CI [{ll95:.3f}; {ul95:.3f}] {test_outcome}")
        plt.title(title.replace('\\n', '\n'))


        plt.plot(es, 0.5, 's', markersize=8, color='black')
        plt.axvline(high_eqbound_d, linestyle='--', color='black')
        plt.axvline(low_eqbound_d, linestyle='--', color='black')
        plt.axvline(0, linestyle='--', color='grey')

        plt.hlines(0.5, ll90, ul90, linewidth=3, color='black')
        plt.hlines(0.5, ll95, ul95, linewidth=1, color='black')

        plt.show()

    if verbose:
        print("TOST results:")
        print(f"Z-value lower bound: {z1:8.2f}\\tp-value lower bound: {p1:.3f}")
        print(f"Z-value upper bound: {z2:8.2f}\\tp-value upper bound: {p2:.3f}")
        print("\\n")
        print("Equivalence bounds (Cohen's d):")
        print(f"low eqbound: {low_eqbound_d:.4f}\\nhigh eqbound: {high_eqbound_d:.4f}")
        print("\\n")
        print("TOST confidence interval:")
        print(f"lower bound {100*(1-alpha*2):.0f}% CI: {ll90:.3f}\\nupper bound {100*(1-alpha*2):.0f}% CI:  {ul90:.3f}")
        print("\\n")
        print("NHST confidence interval:")
        print(f"lower bound {100*(1-alpha):.0f}% CI: {ll95:.3f}\\nupper bound {100*(1-alpha):.0f}% CI:  {ul95:.3f}")
        print("\\n")
        print("Equivalence Test Result:")
        print(f"The equivalence test was {tost_outcome}, Z = {ztost:.3f}, p = {ptost:.3f}, given equivalence bounds of {low_eqbound_d:.3f} and {high_eqbound_d:.3f} and an alpha of {alpha}.")
        print("\\n")
        print("Null Hypothesis Test Result:")
        print(f"The null hypothesis test was {test_outcome}, Z = {z:.3f}, p = {pttest:.3f}, given an alpha of {alpha}.")

        if pttest <= alpha and ptost <= alpha:
            combined_outcome = (f"NHST: reject null significance hypothesis that the effect is equal to 0 \\n"
                                f"TOST: reject null equivalence hypothesis")
        elif pttest < alpha and ptost > alpha:
            combined_outcome = (f"NHST: reject null significance hypothesis that the effect is equal to 0 \\n"
                                f"TOST: don't reject null equivalence hypothesis")
        elif pttest > alpha and ptost <= alpha:
            combined_outcome = (f"NHST: don't reject null significance hypothesis that the effect is equal to 0 \\n"
                                f"TOST: reject null equivalence hypothesis")
        else: # pttest > alpha and ptost > alpha
            combined_outcome = (f"NHST: don't reject null significance hypothesis that the effect is equal to 0 \\n"
                                  f"TOST: don't reject null equivalence hypothesis")

        print("\\n")
        print(combined_outcome)

    results = {
        'ES': es,
        'TOST_Z1': z1,
        'TOST_p1': p1,
        'TOST_Z2': z2,
        'TOST_p2': p2,
        'alpha': alpha,
        'low_eqbound_d': low_eqbound_d,
        'high_eqbound_d': high_eqbound_d,
        'LL_CI_TOST': ll90,
        'UL_CI_TOST': ul90,
        'LL_CI_ZTEST': ll95,
        'UL_CI_ZTEST': ul95,
        'NHST_Z': z,
        'NHST_p': pttest
    }
    return results
