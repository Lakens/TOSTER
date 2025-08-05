import pandas as pd
from scipy import stats
from .anova_summary import anova_summary

def equ_anova(object, eqbound, MET=False, alpha=0.05):
    """
    Equivalence Test for ANOVA Results
    """
    results = anova_summary(object)

    res2 = results.copy()
    # In anova_lm, df is a series of floats, not tuples
    # and the residual df is in the 'Residual' row
    if 'df' in res2.columns and 'F' in res2.columns:
        res2.rename(columns={'df': 'df1', 'F': 'F.value', 'PR(>F)': 'p.null'}, inplace=True)
        res2['df2'] = results.loc['Residual', 'df1']
        res2 = res2.reset_index().rename(columns={'index': 'effect'})
    elif 'Num DF' in res2.columns and 'Den DF' in res2.columns:
        res2.rename(columns={'Num DF': 'df1', 'Den DF': 'df2', 'F Value': 'F.value', 'Pr(>F)': 'p.null'}, inplace=True)
        res2 = res2.reset_index().rename(columns={'index': 'effect'})
    else:
        raise ValueError("Could not determine the columns of the anova table.")

    res2['f2'] = eqbound / (1 - eqbound)
    res2['lambda'] = res2['f2'] * (res2['df1'] + res2['df2'] + 1)

    if MET:
        res2['p.equ'] = 1 - stats.f.cdf(res2['F.value'], res2['df1'], res2['df2'], nc=res2['lambda'])
    else:
        res2['p.equ'] = stats.f.cdf(res2['F.value'], res2['df1'], res2['df2'], nc=res2['lambda'])

    res2['eqbound'] = eqbound

    return res2[['effect', 'df1', 'df2', 'F.value', 'p.null', 'pes', 'eqbound', 'p.equ']]
