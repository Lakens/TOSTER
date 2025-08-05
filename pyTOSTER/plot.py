import numpy as np

def cut_cdf_qi(p, widths=[0.66, 0.95, 1.0]):
    """
    Categorize p-values into confidence levels.
    This is a simplified version of the R function.
    """
    p = np.abs(1 - p * 2)
    # Find which interval each p-value falls into
    cat = np.digitize(p, bins=sorted(widths))
    return np.array(sorted(widths))[cat]
