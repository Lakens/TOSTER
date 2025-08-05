class DataTostOneOptions:
    def __init__(self,
                 vars=None,
                 mu=0,
                 hypothesis="EQU",
                 low_eqbound=-0.5,
                 high_eqbound=0.5,
                 eqbound_type="raw",
                 alpha=0.05,
                 desc=False,
                 plots=False,
                 smd_type="g"):
        self.vars = vars
        self.mu = mu
        self.hypothesis = hypothesis
        self.low_eqbound = low_eqbound
        self.high_eqbound = high_eqbound
        self.eqbound_type = eqbound_type
        self.alpha = alpha
        self.desc = desc
        self.plots = plots
        self.smd_type = smd_type
