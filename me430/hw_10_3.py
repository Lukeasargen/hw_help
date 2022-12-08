import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# Exhaust gas composition

x = 4; y = 10
MW_F = x*12 + y*1
print(f"{MW_F = }")

X_CO2 = 0.120
X_O2 = 0.010
X_H2O = 0.160
X_NO = 200e-6
# Calculate emissions indices
EI_NO = (X_NO/X_CO2)*(x*30/MW_F)
print(f"{EI_NO = }")
