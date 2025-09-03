import numpy as np
from pathlib import Path
import pandas as pd

df_TEs = pd.read_excel('../../Data/TEsummary.xlsx',header=[0],index_col=0)
# Scenarios:
    # 1) Ramp up to 2C, stabilise for 200y @ 2C
        # 1a) Ramp down to 0, stabilise for another 1ky
        # 1b) Remain at 2C, stabilise for another 1ky
    # 2) Ramp up to 4C, stabilise for 200y @ 4C
        # 2a) Ramp down to 0, stabilise for another 1ky
        # 2b) Ramp down to 2C, stabilise for another 1ky
        # 2c) Remain at 4C, stabilise for another 1ky
scenarios = ['1a','1b','2a','2b','2c']
temps = {}
temps['1a'] = np.concatenate((np.arange(0,2.01,step=0.02),
                              np.ones(200)*2,
                              np.arange(0,2,step=0.02)[::-1], 
                              np.zeros(1000)) 
                            )
temps['1b'] = np.concatenate((np.arange(0,2.01,step=0.02),
                         np.ones(1200)*2))

temps['2a'] = np.concatenate((np.arange(0,4.01,step=0.02),
                              np.ones(200)*4,
                              np.arange(0,4,step=0.02)[::-1],
                              np.zeros(1000))
                            )
temps['2b'] = np.concatenate((np.arange(0,4.01,step=0.02),
                              np.ones(200)*4,
                              np.arange(2,4,step=0.02)[::-1],
                              np.ones(1000)*2)
                            )
temps['2c'] = np.concatenate((np.arange(0,4.01,step=0.02),
                         np.ones(1200)*4)
                            )

