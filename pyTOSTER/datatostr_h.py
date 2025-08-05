class DataTostrOptions:
    def __init__(self,
                 pairs=None,
                 cor_type="pearson",
                 hypothesis="EQU",
                 low_eqbound_r=-0.3,
                 high_eqbound_r=0.3,
                 alpha=0.05,
                 desc=False,
                 plots=False):
        self.pairs = pairs
        self.cor_type = cor_type
        self.hypothesis = hypothesis
        self.low_eqbound_r = low_eqbound_r
        self.high_eqbound_r = high_eqbound_r
        self.alpha = alpha
        self.desc = desc
        self.plots = plots
