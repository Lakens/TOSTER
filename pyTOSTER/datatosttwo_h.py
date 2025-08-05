class DataTostTwoOptions:
    def __init__(self,
                 deps=None,
                 group=None,
                 var_equal=False,
                 hypothesis="EQU",
                 low_eqbound=-0.5,
                 high_eqbound=0.5,
                 eqbound_type="raw",
                 alpha=0.05,
                 desc=False,
                 plots=False,
                 descplots=False,
                 smd_type="g"):
        self.deps = deps
        self.group = group
        self.var_equal = var_equal
        self.hypothesis = hypothesis
        self.low_eqbound = low_eqbound
        self.high_eqbound = high_eqbound
        self.eqbound_type = eqbound_type
        self.alpha = alpha
        self.desc = desc
        self.plots = plots
        self.descplots = descplots
        self.smd_type = smd_type
