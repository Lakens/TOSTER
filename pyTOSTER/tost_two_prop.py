import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def tost_two_prop(prop1, prop2, n1, n2, low_eqbound, high_eqbound, alpha=0.05, ci_type="normal", plot=True, verbose=True):
    """
    TOST function for two proportions (raw scores).
    """
    if low_eqbound >= high_eqbound:
        print("Warning: The lower bound is equal to or larger than the upper bound.")
    if n1 < 2 or n2 < 2:
        raise ValueError("The sample size should be larger than 1.")
    if not (0 < alpha < 1):
        raise ValueError("The alpha level should be a positive value between 0 and 1.")
    if not (0 <= prop1 <= 1 and 0 <= prop2 <= 1):
        raise ValueError("The proportion should be a value between 0 and 1.")
    if ci_type not in ["normal", "wilson"]:
        raise ValueError("ci_type must be equal to 'wilson' or 'normal'")

    prop_dif = prop1 - prop2
    prop_se = np.sqrt((prop1 * (1 - prop1)) / n1 + (prop2 * (1 - prop2)) / n2)

    z1 = (prop_dif - low_eqbound) / prop_se
    z2 = (prop_dif - high_eqbound) / prop_se
    z = prop_dif / prop_se

    # p-values for TOST
    p1 = 1 - stats.norm.cdf(z1)
    p2 = stats.norm.cdf(z2)
    ptost = max(p1, p2)
    ztost = z1 if abs(z1) < abs(z2) else z2

    # NHST p-value
    pttest = 2 * (1 - stats.norm.cdf(abs(z)))

    tost_outcome = "significant" if ptost < alpha else "non-significant"
    test_outcome = "significant" if pttest < alpha else "non-significant"

    if ci_type == "normal":
        ll90 = prop_dif - (stats.norm.ppf(1 - alpha) * prop_se)
        ul90 = prop_dif + (stats.norm.ppf(1 - alpha) * prop_se)
        ll95 = prop_dif - (stats.norm.ppf(1 - (alpha / 2)) * prop_se)
        ul95 = prop_dif + (stats.norm.ppf(1 - (alpha / 2)) * prop_se)

    if ci_type == "wilson":
        # This is a direct translation of the R code's Wilson CI,
        # which seems to have a bug. Using proportions instead of p-values.
        # The formula is still likely a simplified version.
        conf_level = 1 - alpha
        delta = prop_dif
        estimate = np.array([prop1, prop2])
        yates = abs(delta) / (1/n1 + 1/n2)

        width90 = stats.norm.ppf(1-alpha) * np.sqrt(np.sum(estimate * (1 - estimate) / np.array([n1, n2]))) + yates * (1/n1 + 1/n2)
        width95 = stats.norm.ppf(1-(alpha/2)) * np.sqrt(np.sum(estimate * (1 - estimate) / np.array([n1, n2]))) + yates * (1/n1 + 1/n2)

        ll90 = max(delta - width90, -1)
        ul90 = min(delta + width90, 1)
        ll95 = max(delta - width95, -1)
        ul95 = min(delta + width95, 1)

    if plot:
        plt.figure()
        plt.ylim(0, 1)
        plot_min = min(ll95, low_eqbound) - max(ul95 - ll95, high_eqbound - low_eqbound) / 10
        plot_max = max(ul95, high_eqbound) + max(ul95 - ll95, high_eqbound - low_eqbound) / 10
        plt.xlim(plot_min, plot_max)

        plt.yticks([])
        plt.xlabel("Proportion Difference")
        title = (f"Equivalence bounds {low_eqbound:.3f} and {high_eqbound:.3f}\\n"
                 f"Proportion Difference = {prop_dif:.3f}\\n"
                 f"TOST: {100*(1-alpha*2):.0f}% CI [{ll90:.3f}; {ul90:.3f}] {tost_outcome}\\n"
                 f"NHST: {100*(1-alpha):.0f}% CI [{ll95:.3f}; {ul95:.3f}] {test_outcome}")
        plt.title(title.replace('\\n', '\n'))

        plt.plot(prop_dif, 0.5, 's', markersize=8, color='black')
        plt.axvline(high_eqbound, linestyle='--', color='black')
        plt.axvline(low_eqbound, linestyle='--', color='black')
        plt.axvline(0, linestyle='--', color='grey')

        plt.hlines(0.5, ll90, ul90, linewidth=3, color='black')
        plt.hlines(0.5, ll95, ul95, linewidth=1, color='black')

        plt.show()

    if verbose:
        print("TOST results:")
        print(f"Z-value lower bound: {z1:8.2f}\\tp-value lower bound: {p1:.3f}")
        print(f"Z-value upper bound: {z2:8.2f}\\tp-value upper bound: {p2:.3f}")
        print("\\n")
        print("Equivalence bounds:")
        print(f"low eqbound: {low_eqbound:.4f}\\nhigh eqbound: {high_eqbound:.4f}")
        print("\\n")
        print("TOST confidence interval:")
        print(f"lower bound {100*(1-alpha*2):.0f}% CI: {ll90:.3f}\\nupper bound {100*(1-alpha*2):.0f}% CI:  {ul90:.3f}")
        print("\\n")
        print("NHST confidence interval:")
        print(f"lower bound {100*(1-alpha):.0f}% CI: {ll95:.3f}\\nupper bound {100*(1-alpha):.0f}% CI:  {ul95:.3f}")
        print("\\n")
        print("Equivalence Test based on Fisher's exact z-test Result:")
        print(f"The equivalence test was {tost_outcome}, Z = {ztost:.3f}, p = {ptost:.3f}, given equivalence bounds of {low_eqbound:.3f} and {high_eqbound:.3f} and an alpha of {alpha}.")
        print("\\n")
        print("Null-Hypothesis Fisher's exact z-test Result:")
        print(f"The null hypothesis test was {test_outcome}, Z = {z:.3f}, p = {pttest:.3f}, given an alpha of {alpha}.")

        if pttest <= alpha and ptost <= alpha:
            combined_outcome = ("NHST: reject null significance hypothesis\\n"
                                "TOST: reject null equivalence hypothesis")
        elif pttest < alpha and ptost > alpha:
            combined_outcome = ("NHST: reject null significance hypothesis\\n"
                                "TOST: don't reject null equivalence hypothesis")
        elif pttest > alpha and ptost <= alpha:
            combined_outcome = ("NHST: don't reject null significance hypothesis that the effect is equal to 0\\n"
                                "TOST: reject null equivalence hypothesis")
        else: # pttest > alpha and ptost > alpha
            combined_outcome = ("NHST: don't reject null significance hypothesis that the effect is equal to 0\\n"
                                  "TOST: don't reject null equivalence hypothesis")

        print("\\n")
        print(combined_outcome)

    results = {
        'dif': prop_dif,
        'TOST_z1': z1,
        'TOST_p1': p1,
        'TOST_z2': z2,
        'TOST_p2': p2,
        'alpha': alpha,
        'low_eqbound': low_eqbound,
        'high_eqbound': high_eqbound,
        'LL_CI_TOST': ll90,
        'UL_CI_TOST': ul90,
        'LL_CI_ZTEST': ll95,
        'UL_CI_ZTEST': ul95,
        'NHST_z': z,
        'NHST_p': pttest
    }
    return results
