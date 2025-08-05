import numpy as np
from .two_proportions import twoprop_test

class DataTostTwoProp:
    def __init__(self, data, var, group, level, alpha=0.05,
                 low_eqbound=None, high_eqbound=None,
                 hypothesis='EQU'):
        self.data = data
        self.var = var
        self.group = group
        self.level = level
        self.alpha = alpha
        self.low_eqbound = low_eqbound
        self.high_eqbound = high_eqbound
        self.hypothesis = hypothesis
        self.results = {}

    def run(self):
        dep_data = self.data[self.var]
        group_data = self.data[self.group]

        levels = np.unique(group_data)
        if len(levels) != 2:
            raise ValueError("Grouping variable must have exactly 2 levels.")

        counts1 = np.sum(dep_data[group_data == levels[0]] == self.level)
        n1 = len(dep_data[group_data == levels[0]])
        prop1 = counts1 / n1

        counts2 = np.sum(dep_data[group_data == levels[1]] == self.level)
        n2 = len(dep_data[group_data == levels[1]])
        prop2 = counts2 / n2

        if self.hypothesis == 'EQU':
            alternative = 'equivalence'
        elif self.hypothesis == 'MET':
            alternative = 'minimal.effect'
        else:
            alternative = 'two.sided'

        tost_res = twoprop_test(p1=prop1, p2=prop2,
                                n1=n1, n2=n2,
                                alternative=alternative,
                                null=[self.low_eqbound, self.high_eqbound],
                                alpha=self.alpha)

        self.results = tost_res
        return self.results
