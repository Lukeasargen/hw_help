import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

def rho_ratio_from_P_ratio(P_ratio, gamma=1.4):  # rho2/rho1
    k = (gamma+1)/(gamma-1)
    return (1 + k*P_ratio)/(k + P_ratio)

# Constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

P1 = 101325  # Pa
T1 = 298.15  # K
Vs = 880  # m/s, shock speed
P_ratio = 7.375  # P2/P1

a1 = np.sqrt(gamma*R*T1)
rho_ratio = rho_ratio_from_P_ratio(P_ratio, gamma=gamma)
T_ratio = P_ratio/rho_ratio  # T2/T1
Vas = Vs*(1 - 1/rho_ratio)
Ms = np.sqrt(1 + ((gamma+1)/(2*gamma))*(P_ratio-1))
Vs = Ms*a1

print(f"{rho_ratio = }")
print(f"{T_ratio = }")
print(f"{Vas = } m/s")
print(f"{Ms = }")
print(f"{Vs = } m/s")
