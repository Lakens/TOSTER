class DataTostTwoPropOptions:
    def __init__(self,
                 var=None,
                 level=None,
                 group=None,
                 hypothesis="EQU",
                 low_eqbound=-0.1,
                 high_eqbound=0.1,
                 alpha=0.05,
                 desc=False,
                 plot=False):
        self.var = var
        self.level = level
        self.group = group
        self.hypothesis = hypothesis
        self.low_eqbound = low_eqbound
        self.high_eqbound = high_eqbound
        self.alpha = alpha
        self.desc = desc
        self.plot = plot
