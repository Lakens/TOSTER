def as_htest(tost_res):
    """
    Converts a TOSTER result object to a dictionary similar to R's htest object.
    """
    if not isinstance(tost_res, dict) or 'TOST' not in tost_res:
        raise TypeError("Input must be a TOSTER result dictionary.")

    tost_df = tost_res.get('TOST')

    if tost_df is None:
        # Handle cases where TOST table is not present, e.g., from boot_t_test
        return {
            'statistic': tost_res.get('statistic'),
            'p.value': tost_res.get('p.value'),
            'conf.int': tost_res.get('conf.int'),
            'estimate': tost_res.get('estimate'),
            'null.value': tost_res.get('null.value'),
            'alternative': tost_res.get('alternative'),
            'method': tost_res.get('method')
        }

    if 'minimal.effect' in tost_res.get('alternative', ''):
        tostp = tost_df.loc[['TOST Lower', 'TOST Upper']].sort_values('p.value').iloc[0]
    else: # equivalence
        tostp = tost_df.loc[['TOST Lower', 'TOST Upper']].sort_values('p.value', ascending=False).iloc[0]

    htest = {
        'statistic': tostp['t'],
        'parameter': tostp['df'],
        'p.value': tostp['p.value'],
        'estimate': tost_res['effsize'].loc['Raw', 'estimate'],
        'null.value': [tost_res['eqb'].loc['Raw', 'low_eq'], tost_res['eqb'].loc['Raw', 'high_eq']],
        'alternative': tost_res.get('alternative', 'two.sided'),
        'method': tost_res.get('method', ''),
        'conf.int': [tost_res['effsize'].loc['Raw', 'lower.ci'], tost_res['effsize'].loc['Raw', 'upper.ci']]
    }

    return htest
