import pandas as pd

def anova_summary(anova_table):
    """
    Calculates partial eta squared for an ANOVA table.

    Parameters
    ----------
    anova_table : pandas.DataFrame
        An ANOVA table from statsmodels.

    Returns
    -------
    pandas.DataFrame
        The ANOVA table with a partial eta squared column.
    """
    if not isinstance(anova_table, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame from statsmodels ANOVA.")

    # anova_lm has df, F, and PR(>F)
    # AnovaRM has F Value, Num DF, Den DF, Pr(>F)

    # Make a copy to avoid SettingWithCopyWarning
    anova_table = anova_table.copy()

    if 'F' in anova_table.columns and 'df' in anova_table.columns:
        # From stats.anova_lm
        df_num = anova_table['df']
        df_den = anova_table.loc['Residual', 'df']
        pes = (anova_table['F'] * df_num) / (anova_table['F'] * df_num + df_den)
        anova_table['pes'] = pes
    elif 'F Value' in anova_table.columns and 'Num DF' in anova_table.columns and 'Den DF' in anova_table.columns:
        # From AnovaRM
        pes = (anova_table['F Value'] * anova_table['Num DF']) / (anova_table['F Value'] * anova_table['Num DF'] + anova_table['Den DF'])
        anova_table['pes'] = pes
    elif 'F.value' in anova_table.columns and 'df1' in anova_table.columns and 'df2' in anova_table.columns:
        # From other sources, like the output of the other functions in this package
        pes = (anova_table['F.value'] * anova_table['df1']) / (anova_table['F.value'] * anova_table['df1'] + anova_table['df2'])
        anova_table['pes'] = pes
    else:
        print("Warning: Could not determine F, df1, and df2 columns to calculate partial eta squared.")

    return anova_table
