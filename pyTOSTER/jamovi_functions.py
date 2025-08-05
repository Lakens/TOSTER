import numpy as np
import matplotlib.pyplot as plt

def plot_tost_jam(x, type="cd", estimates=["raw", "SMD"], ggtheme=None):

    smd_label = x['smd']['smd_label']
    low_eqd = x['eqb'].loc[x['eqb']['type'] == smd_label, 'low_eq'].iloc[0]
    high_eqd = x['eqb'].loc[x['eqb']['type'] == smd_label, 'high_eq'].iloc[0]

    low_eqt = x['eqb'].loc[x['eqb']['type'] == 'Raw', 'low_eq'].iloc[0]
    high_eqt = x['eqb'].loc[x['eqb']['type'] == 'Raw', 'high_eq'].iloc[0]


    fig, axes = plt.subplots(len(estimates), 1, figsize=(6, 4*len(estimates)))
    if len(estimates) == 1:
        axes = [axes]

    # SMD plot
    if "SMD" in estimates:
        ax = axes[0]
        smd_est = x['smd']['d']
        smd_ci = [x['smd']['dlow'], x['smd']['dhigh']]
        ax.plot(smd_est, 0, 'o', markersize=8)
        ax.plot(smd_ci, [0, 0], 'k-', linewidth=2)
        ax.axvline(low_eqd, linestyle='--', color='red')
        ax.axvline(high_eqd, linestyle='--', color='red')
        ax.set_title(smd_label)
        ax.set_yticks([])

    # Raw plot
    if "raw" in estimates:
        ax = axes[len(estimates)-1]
        raw_est = x['effsize'].loc['Raw', 'estimate']
        raw_ci = [x['effsize'].loc['Raw', 'lower.ci'], x['effsize'].loc['Raw', 'upper.ci']]
        ax.plot(raw_est, 0, 'o', markersize=8)
        ax.plot(raw_ci, [0, 0], 'k-', linewidth=2)
        ax.axvline(low_eqt, linestyle='--', color='red')
        ax.axvline(high_eqt, linestyle='--', color='red')
        ax.set_title("Raw")
        ax.set_yticks([])

    plt.tight_layout()
    plt.show()
