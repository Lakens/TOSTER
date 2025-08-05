import pandas as pd
from .htest_helpers import as_htest, describe_htest
from .jamovi_functions import plot_tost_jam
import warnings

def print_tostt(x, digits=4):
    effsize = x['effsize']
    tost = x['TOST']

    tost['p.value'] = [f"< 0.001" if p < 0.001 else round(p, 3) for p in tost['p.value']]
    effsize['CI'] = [f"[{round(l, digits)}, {round(u, digits)}]" for l, u in zip(effsize['lower.ci'], effsize['upper.ci'])]

    print(x['method'])
    print(x['decision']['TOST'])
    print(x['decision']['ttest'])
    print(x['decision']['combined'])
    print("\nTOST Results")
    print(tost)
    print("\nEffect Sizes")
    print(effsize[['estimate', 'SE', 'CI', 'conf.level']])

def describe_tostt(x, digits=3):
    return describe_htest(as_htest(x), digits=digits)

def plot_tostt(x, type="simple", estimates=["raw", "SMD"]):
    if type == "simple":
        plot_tost_jam(x, estimates=estimates)
    else:
        warnings.warn("Only 'simple' plot type is supported in this Python version.")
        plot_tost_jam(x, estimates=estimates)
