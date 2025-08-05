class DataTostPairedOptions:
    def __init__(self,
                 pair1=None,
                 pair2=None,
                 hypothesis="EQU",
                 low_eqbound=-0.5,
                 high_eqbound=0.5,
                 eqbound_type="raw",
                 alpha=0.05,
                 desc=False,
                 plots=False,
                 smd_type="g"):
        self.pair1 = pair1
        self.pair2 = pair2
        self.hypothesis = hypothesis
        self.low_eqbound = low_eqbound
        self.high_eqbound = high_eqbound
        self.eqbound_type = eqbound_type
        self.alpha = alpha
        self.desc = desc
        self.plots = plots
        self.smd_type = smd_type
