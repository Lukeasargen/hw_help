import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# Constants for air
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K



# premixed turbulent isoctane-air
equiv_ratio = 0.9
T_u = 500  # K, reactant temperature
P_u = 6*101325  # Pa, pressure
u_rms = 8  # m/s, rms turbulence velocity
l_0 = 6*0.001  # m, turbulence integral length scale
dynamic_viscosity = 270.1e-7  # Ns/m^2, Table C.1 at 500K
alpha = (1/6)*56.7e-6  # m^2/s, Table C.1 at 500K, thermal diffusivity
# Metghalchi & Keck flame-speed correlation
T_ref = 300  # K
P_ref = 1*101325  # Pa
phi_m = 1.15
B_m = 26.32 /100  # m/s
B_2 = -84.72 /100  # m/s
Y_dil = 0  # no dilution

# 1a
tau_0 = l_0/u_rms
print(f"{tau_0 = } s")

# 1b
rho_air = P_u/(R*T_u)
print(f"{rho_air = } kg/m^3")
Re_l_0 = (rho_air*u_rms*l_0)/dynamic_viscosity
print(f"{Re_l_0 = }")

# 1c
S_L_ref = B_m + B_2*(equiv_ratio-phi_m)**2
print(f"{S_L_ref = } m/s")
gamma_exponent = 2.18 - 0.8*(equiv_ratio-1)
print(f"{gamma_exponent = }")
beta_exponent = -0.16 + 0.22*(equiv_ratio-1)
print(f"{beta_exponent = }")
S_L = S_L_ref*((T_u/T_ref)**gamma_exponent)*((P_u/P_ref)**beta_exponent)*(1-2.1*Y_dil)
thickness = 2*alpha/(S_L)  # m
print(f"{S_L = } m/s")
print(f"{thickness = } m")
tau_chem = thickness/S_L
print(f"{tau_chem = } s")

# 1d
Da = tau_0/tau_chem
print(f"{Da = }")
# Da = (l_0/thickness)*(S_L/u_rms)
# print(f"{Da = }")

# 1e
print(f"{u_rms/S_L = }")
print(f"{thickness/l_0 = }")
