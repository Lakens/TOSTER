import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def gg_curv_t(data_list,
              levels=[0.68, 0.90, 0.95, 0.999],
              xaxis='Range of Values',
              yaxis1='two-tailed p-value',
              yaxis2='Confidence Interval (%)',
              color='black',
              fill='skyblue',
              alpha_shade=0.5):

    data = data_list[0]

    # Sort levels for shading
    ci_shade1 = sorted(levels, reverse=True)

    # Create intervals for points
    interval_data = []
    for level in ci_shade1:
        idx = np.abs(data['intrvl.level'] - level).idxmin()
        interval_data.append({'levels': level, 'limits': data.loc[idx, 'lower.limit']})
        interval_data.append({'levels': level, 'limits': data.loc[idx, 'upper.limit']})

    interval_df = pd.DataFrame(interval_data)

    fig, ax1 = plt.subplots()

    # Plot consonance curve
    ax1.plot(data['lower.limit'], data['pvalue'], color=color)
    ax1.plot(data['upper.limit'], data['pvalue'], color=color)
    ax1.fill_between(data['lower.limit'], 0, data['pvalue'], color=fill, alpha=alpha_shade)
    ax1.fill_between(data['upper.limit'], 0, data['pvalue'], color=fill, alpha=alpha_shade)

    # Plot interval points
    ax1.scatter(interval_df['limits'], 1 - interval_df['levels'], color=color, marker='d')

    # Plot interval lines
    for level in ci_shade1:
        subset = interval_df[interval_df['levels'] == level]
        ax1.plot(subset['limits'], 1 - subset['levels'], color=color, linewidth=0.5)

    ax1.set_xlabel(xaxis)
    ax1.set_ylabel(yaxis1)

    # Secondary y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel(yaxis2)
    ax2.set_ylim(0, 100)
    ax1.set_ylim(0, 1)

    # Set ticks for secondary axis
    ax2.set_yticks(np.arange(0, 101, 10))
    ax1.set_yticks(np.arange(0, 1.1, 0.1))

    plt.show()
