import numpy as np
from .cor_test import cor_test

class DataTostr:
    def __init__(self, data, pairs, alpha=0.05,
                 low_eqbound_r=None, high_eqbound_r=None,
                 cor_type='pearson', hypothesis='EQU'):
        self.data = data
        self.pairs = pairs
        self.alpha = alpha
        self.low_eqbound_r = low_eqbound_r
        self.high_eqbound_r = high_eqbound_r
        self.cor_type = cor_type
        self.hypothesis = hypothesis
        self.results = {}

    def run(self):
        for pair in self.pairs:
            pair1_data = self.data[pair[0]]
            pair2_data = self.data[pair[1]]

            mask = ~np.isnan(pair1_data) & ~np.isnan(pair2_data)
            pair1_data = pair1_data[mask]
            pair2_data = pair2_data[mask]

            if self.hypothesis == 'EQU':
                alternative = 'equivalence'
            elif self.hypothesis == 'MET':
                alternative = 'minimal.effect'
            else:
                alternative = 'two.sided'

            res = cor_test(x=pair1_data, y=pair2_data,
                           method=self.cor_type,
                           alternative=alternative,
                           null=[self.low_eqbound_r, self.high_eqbound_r],
                           alpha=self.alpha)

            self.results[f"{pair[0]}_{pair[1]}"] = res

        return self.results
