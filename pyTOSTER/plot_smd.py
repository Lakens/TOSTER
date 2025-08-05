import numpy as np
import pandas as pd
from .gg_curv_t import gg_curv_t
from .plot_scripts import d_curv_raw

def plot_smd(d, df, lambda_val=None, sigma=None, smd_ci="t",
             smd_label="SMD", type=["c", "cd"],
             levels=[0.5, 0.9, 0.95, 0.999]):

    if smd_ci == "nct":
        raise NotImplementedError("nct method not supported for this function at this time.")

    if smd_ci == "goulet" and lambda_val is None:
        raise ValueError("lambda_val must be provided when smd_ci is 'goulet'")

    if smd_ci in ["t", "z"] and sigma is None:
        raise ValueError("sigma must be provided when smd_ci is 't' or 'z'")

    dat = d_curv_raw(d=d, df=df, ncp=lambda_val, sigma=sigma, smd_ci=smd_ci)

    return gg_curv_t(dat, type=type, levels=levels,
                     xaxis=smd_label,
                     yaxis1="p-value",
                     yaxis2="Confidence Interval (%)")
