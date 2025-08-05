import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def tost_one(m, mu, sd, n, low_eqbound_d, high_eqbound_d, alpha=0.05, plot=True, verbose=True):
    """
    TOST function for a one-sample t-test (Cohen's d).

    Parameters
    ----------
    m : float
        mean.
    mu : float
        value to compare against.
    sd : float
        standard deviation.
    n : int
        sample size.
    low_eqbound_d : float
        lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's d).
    high_eqbound_d : float
        upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's d).
    alpha : float, optional
        alpha level (default = 0.05).
    plot : bool, optional
        set whether results should be plotted (default is True).
    verbose : bool, optional
        logical variable indicating whether text output should be generated (default is True).

    Returns
    -------
    dict
        A dictionary containing the TOST results.
    """

    if low_eqbound_d >= high_eqbound_d:
        print("Warning: The lower bound is equal to or larger than the upper bound.")
    if n < 2:
        raise ValueError("The sample size should be larger than 1.")
    if not (0 < alpha < 1):
        raise ValueError("The alpha level should be a positive value between 0 and 1.")
    if sd <= 0:
        raise ValueError("The standard deviation should be a positive value.")

    low_eqbound = low_eqbound_d * sd
    high_eqbound = high_eqbound_d * sd

    return tost_one_raw(m, mu, sd, n, low_eqbound, high_eqbound, alpha, plot, verbose,
                        low_eqbound_d=low_eqbound_d, high_eqbound_d=high_eqbound_d)


def tost_one_raw(m, mu, sd, n, low_eqbound, high_eqbound, alpha=0.05, plot=True, verbose=True, low_eqbound_d=None, high_eqbound_d=None):
    """
    TOST function for a one-sample t-test (raw scores).
    """

    if low_eqbound >= high_eqbound:
        print("Warning: The lower bound is equal to or larger than the upper bound.")
    if n < 2:
        raise ValueError("The sample size should be larger than 1.")
    if not (0 < alpha < 1):
        raise ValueError("The alpha level should be a positive value between 0 and 1.")
    if sd <= 0:
        raise ValueError("The standard deviation should be a positive value.")

    degree_f = n - 1

    t1 = (m - mu - low_eqbound) / (sd / np.sqrt(n))
    p1 = 1 - stats.t.cdf(t1, df=degree_f)
    t2 = (m - mu - high_eqbound) / (sd / np.sqrt(n))
    p2 = stats.t.cdf(t2, df=degree_f)

    t = (m - mu) / (sd / np.sqrt(n))
    pttest = 2 * (1 - stats.t.cdf(abs(t), df=degree_f))

    ll90 = (m - mu) - stats.t.ppf(1 - alpha, df=degree_f) * (sd / np.sqrt(n))
    ul90 = (m - mu) + stats.t.ppf(1 - alpha, df=degree_f) * (sd / np.sqrt(n))
    ll95 = (m - mu) - stats.t.ppf(1 - (alpha / 2), df=degree_f) * (sd / np.sqrt(n))
    ul95 = (m - mu) + stats.t.ppf(1 - (alpha / 2), df=degree_f) * (sd / np.sqrt(n))

    ptost = max(p1, p2)
    ttost = t1 if abs(t1) < abs(t2) else t2
    diff = m - mu

    test_outcome = "significant" if pttest < alpha else "non-significant"
    tost_outcome = "significant" if ptost < alpha else "non-significant"

    if plot:
        plt.figure()
        plt.ylim(0, 1)
        plot_min = min(ll95, low_eqbound, diff) - max(ul95 - ll95, high_eqbound - low_eqbound, diff) / 10
        plot_max = max(ul95, high_eqbound, diff) + max(ul95 - ll95, high_eqbound - low_eqbound, diff) / 10
        plt.xlim(plot_min, plot_max)

        plt.yticks([])
        plt.xlabel("Mean Difference")
        title = (f"Equivalence bounds {low_eqbound:.3f} and {high_eqbound:.3f}\\n"
                 f"Mean difference = {diff:.3f}\\n"
                 f"TOST: {100*(1-alpha*2):.0f}% CI [{ll90:.3f}; {ul90:.3f}] {tost_outcome}\\n"
                 f"NHST: {100*(1-alpha):.0f}% CI [{ll95:.3f}; {ul95:.3f}] {test_outcome}")
        plt.title(title.replace('\\n', '\n'))

        plt.plot(diff, 0.5, 's', markersize=8, color='black')
        plt.axvline(high_eqbound, linestyle='--', color='black')
        plt.axvline(low_eqbound, linestyle='--', color='black')
        plt.axvline(0, linestyle='--', color='grey')

        plt.hlines(0.5, ll90, ul90, linewidth=3, color='black')
        plt.hlines(0.5, ll95, ul95, linewidth=1, color='black')

        plt.show()

    if verbose:
        print("TOST results:")
        print(f"t-value lower bound: {t1:8.2f}\\tp-value lower bound: {p1:.3f}")
        print(f"t-value upper bound: {t2:8.2f}\\tp-value upper bound: {p2:.3f}")
        print(f"degrees of freedom : {degree_f:.2f}")
        print("\\n")
        if low_eqbound_d is not None and high_eqbound_d is not None:
            print("Equivalence bounds (Cohen's d):")
            print(f"low eqbound: {low_eqbound_d:.4f}\\nhigh eqbound: {high_eqbound_d:.4f}")
            print("\\n")
        print("Equivalence bounds (raw scores):")
        print(f"low eqbound: {low_eqbound:.4f}\\nhigh eqbound: {high_eqbound:.4f}")
        print("\\n")
        print("TOST confidence interval:")
        print(f"lower bound {100*(1-alpha*2):.0f}% CI: {ll90:.3f}\\nupper bound {100*(1-alpha*2):.0f}% CI:  {ul90:.3f}")
        print("\\n")
        print("NHST confidence interval:")
        print(f"lower bound {100*(1-alpha):.0f}% CI: {ll95:.3f}\\nupper bound {100*(1-alpha):.0f}% CI:  {ul95:.3f}")
        print("\\n")
        print("Equivalence Test Result:")
        print(f"The equivalence test was {tost_outcome}, t({degree_f:.2f}) = {ttost:.3f}, p = {ptost:.3f}, given equivalence bounds of {low_eqbound:.3f} and {high_eqbound:.3f} (on a raw scale) and an alpha of {alpha}.")
        print("\\n")
        print("Null Hypothesis Test Result:")
        print(f"The null hypothesis test was {test_outcome}, t({degree_f:.2f}) = {t:.3f}, p = {pttest:.3f}, given an alpha of {alpha}.")

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
        'diff': diff,
        'TOST_t1': t1,
        'TOST_p1': p1,
        'TOST_t2': t2,
        'TOST_p2': p2,
        'TOST_df': degree_f,
        'alpha': alpha,
        'low_eqbound': low_eqbound,
        'high_eqbound': high_eqbound,
        'LL_CI_TOST': ll90,
        'UL_CI_TOST': ul90,
        'LL_CI_TTEST': ll95,
        'UL_CI_TTEST': ul95,
        'NHST_t': t,
        'NHST_p': pttest
    }
    if low_eqbound_d is not None:
        results['low_eqbound_d'] = low_eqbound_d
    if high_eqbound_d is not None:
        results['high_eqbound_d'] = high_eqbound_d

    return results
