import pandas as pd
import numpy as np

class new:
    def __init__(self, N):
        empty = np.zeros(N)
        my_results = {
            "M_re": empty,
            "M_x": empty,
            "M_ch": empty,
            "MgO": empty,
            "Th": empty,
            "Al2O3": empty,
            "Fe3+": empty,
            "H2O": empty,
            "S_m": empty,
            "S_f": empty,
            "fluid": empty,
        }
        self.df = pd.DataFrame(data=my_results)
