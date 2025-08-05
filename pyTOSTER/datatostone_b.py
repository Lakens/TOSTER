import numpy as np
from .t_tost import t_tost

class DataTostOne:
    def __init__(self, data, vars, alpha=0.05, mu=0,
                 low_eqbound=None, high_eqbound=None,
                 eqbound_type='raw', hypothesis='EQU', smd_type='g'):
        self.data = data
        self.vars = vars
        self.alpha = alpha
        self.mu = mu
        self.low_eqbound = low_eqbound
        self.high_eqbound = high_eqbound
        self.eqbound_type = eqbound_type
        self.hypothesis = hypothesis
        self.smd_type = smd_type
        self.results = {}

    def run(self):
        for var_name in self.vars:
            var_data = self.data[var_name]
            var_data = var_data[~np.isnan(var_data)]

            bias_c = True if self.smd_type == 'g' else False

            tost_res = t_tost(x=var_data,
                               hypothesis=self.hypothesis,
                               low_eqbound=self.low_eqbound,
                               high_eqbound=self.high_eqbound,
                               eqbound_type=self.eqbound_type,
                               alpha=self.alpha,
                               bias_correction=bias_c,
                               mu=self.mu)

            self.results[var_name] = tost_res

        return self.results
