import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# Exhaust gas composition
x = 8; y = 18
MW_F = x*12 + y*1
print(f"{MW_F = }")
# All water removed prior to measurements
X_CO2 = 0.0844
X_O2 = 0.0879
X_CO = 44e-6
X_NO = 76e-6
X_UHC = 15e-6
# Calculate emissions indices
print(X_CO2+X_CO+X_UHC)
EI_CO = (X_CO/(X_CO2+X_CO+X_UHC))*(x*28/MW_F)
print(f"{EI_CO = }")
EI_NO = (X_NO/(X_CO2+X_CO+X_UHC))*(x*30/MW_F)
print(f"{EI_NO = }")
EI_UHC = (X_UHC/(X_CO2+X_CO+X_UHC))*(x*16/MW_F)
print(f"{EI_UHC = }")
