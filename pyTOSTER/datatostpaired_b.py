import numpy as np
from .t_tost import t_tost

class DataTostPaired:
    def __init__(self, data, pair1, pair2, alpha=0.05,
                 low_eqbound=None, high_eqbound=None,
                 eqbound_type='raw', hypothesis='EQU', smd_type='g'):
        self.data = data
        self.pair1 = pair1
        self.pair2 = pair2
        self.alpha = alpha
        self.low_eqbound = low_eqbound
        self.high_eqbound = high_eqbound
        self.eqbound_type = eqbound_type
        self.hypothesis = hypothesis
        self.smd_type = smd_type
        self.results = {}

    def run(self):
        pair1_data = self.data[self.pair1]
        pair2_data = self.data[self.pair2]

        mask = ~np.isnan(pair1_data) & ~np.isnan(pair2_data)
        pair1_data = pair1_data[mask]
        pair2_data = pair2_data[mask]

        bias_c = True if self.smd_type == 'g' else False

        tost_res = t_tost(x=pair1_data, y=pair2_data,
                           paired=True,
                           hypothesis=self.hypothesis,
                           low_eqbound=self.low_eqbound,
                           high_eqbound=self.high_eqbound,
                           eqbound_type=self.eqbound_type,
                           alpha=self.alpha,
                           bias_correction=bias_c)

        self.results = tost_res
        return self.results
