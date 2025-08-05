import numpy as np
from .t_tost import t_tost

class DataTostTwo:
    def __init__(self, data, deps, group, alpha=0.05,
                 low_eqbound=None, high_eqbound=None,
                 eqbound_type='raw', hypothesis='EQU', smd_type='g',
                 var_equal=False):
        self.data = data
        self.deps = deps
        self.group = group
        self.alpha = alpha
        self.low_eqbound = low_eqbound
        self.high_eqbound = high_eqbound
        self.eqbound_type = eqbound_type
        self.hypothesis = hypothesis
        self.smd_type = smd_type
        self.var_equal = var_equal
        self.results = {}

    def run(self):
        for dep_name in self.deps:
            dep_data = self.data[dep_name]
            group_data = self.data[self.group]

            levels = np.unique(group_data)
            if len(levels) != 2:
                raise ValueError("Grouping variable must have exactly 2 levels.")

            data1 = dep_data[group_data == levels[0]]
            data2 = dep_data[group_data == levels[1]]

            bias_c = True if self.smd_type == 'g' else False

            tost_res = t_tost(x=data1, y=data2,
                               paired=False,
                               var_equal=self.var_equal,
                               hypothesis=self.hypothesis,
                               low_eqbound=self.low_eqbound,
                               high_eqbound=self.high_eqbound,
                               eqbound_type=self.eqbound_type,
                               alpha=self.alpha,
                               bias_correction=bias_c)

            self.results[dep_name] = tost_res

        return self.results
