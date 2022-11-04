import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# Metghalchi & Keck flame-speed correlation
phi_m = 1.05
B_m = 50  # cm/s
B_2 = -100  # cm/s
Y_dil = 0  # no dilution

alpha_ref = 20e-6  # m^2/s, thermal diffusivity
T_ref = 300  # K
P_ref = 1  # atm

# 1a) reference flame speed vs equiv ratio
equiv_ratio = np.arange(0.5, 1.65, 0.05)
S_L_ref = B_m + B_2*(equiv_ratio-phi_m)**2
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
# ax.plot(equiv_ratios, S_L_ref)
# ax.set_title("Reference laminar flame speed vs equivalence ratio.")
# ax.set_ylabel("Reference Laminar Flame Speed [cm/s]")
# ax.set_xlabel("Equivalence Ratio")
# ax.grid()
# fig.tight_layout()

# 2b) flame speed vs equiv ratio at conditions
P_u = 1  # atm
T_u = 400  # K
gamma_exponent = 2.18 - 0.8*(equiv_ratio-1)
beta_exponent = -0.16 + 0.22*(equiv_ratio-1)
# flame speed, eq 8.33, cm/s
S_L = S_L_ref*((T_u/T_ref)**gamma_exponent)*((P_u/P_ref)**beta_exponent)*(1-2.1*Y_dil)
# estimate thermal diffusivity at conditions
alpha = alpha_ref*((T_u/T_ref)**1.6)*(P_ref/P_u)
# flame thickness, 10mm = cm
thickness = 1000 * 2*alpha/(S_L/10)  # mm
# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(6,3))
# ax1.plot(equiv_ratio, S_L)
# ax1.set_title("Laminar flame speed vs equivalence ratio.")
# ax1.set_ylabel("Laminar Flame Speed [cm/s]")
# ax1.set_xlabel("Equivalence Ratio")
# ax1.grid()
# ax2.plot(equiv_ratio, thickness)
# ax2.set_title("Laminar flame thickness vs equivalence ratio.")
# ax2.set_ylabel("Laminar Flame Thickness [mm]")
# ax2.set_xlabel("Equivalence Ratio")
# ax2.grid()
# fig.tight_layout()

# 1c) laminar flame speed vs temperature
T_u = np.arange(300, 1300, 100)
equiv_ratio = 1
P_u = 1
S_L_ref = B_m + B_2*(equiv_ratio-phi_m)**2
gamma_exponent = 2.18 - 0.8*(equiv_ratio-1)
beta_exponent = -0.16 + 0.22*(equiv_ratio-1)
S_L = S_L_ref*((T_u/T_ref)**gamma_exponent)*((P_u/P_ref)**beta_exponent)*(1-2.1*Y_dil)
alpha = alpha_ref*((T_u/T_ref)**1.6)*(P_ref/P_u)
thickness = 1000 * 2*alpha/(S_L/10)  # mm
# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(7,3))
# ax1.plot(T_u, S_L)
# ax1.set_title("Laminar flame speed vs unburned gas temperature.")
# ax1.set_ylabel("Laminar Flame Speed [cm/s]")
# ax1.set_xlabel("Unburned Gas Temperature [K]")
# ax1.grid()
# ax2.plot(T_u, thickness)
# ax2.set_title("Laminar flame thickness vs unburned gas temperature.")
# ax2.set_ylabel("Laminar Flame Thickness [mm]")
# ax2.set_xlabel("Unburned Gas Temperature [K]")
# ax2.grid()
# fig.tight_layout()

# 1d) laminar flame speed vs pressure
P_u = np.array([0.5, 1, 2, 4, 7, 16, 32, 64])
T_u = 400
equiv_ratio = 1
S_L_ref = B_m + B_2*(equiv_ratio-phi_m)**2
gamma_exponent = 2.18 - 0.8*(equiv_ratio-1)
beta_exponent = -0.16 + 0.22*(equiv_ratio-1)
S_L = S_L_ref*((T_u/T_ref)**gamma_exponent)*((P_u/P_ref)**beta_exponent)*(1-2.1*Y_dil)
alpha = alpha_ref*((T_u/T_ref)**1.6)*(P_ref/P_u)
thickness = 1000 * 2*alpha/(S_L/10)  # mm
# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(7,3))
# ax1.plot(P_u, S_L)
# ax1.set_title("Laminar flame speed vs pressure.")
# ax1.set_ylabel("Laminar Flame Speed [cm/s]")
# ax1.set_xlabel("Pressure [atm]")
# ax1.set_xscale('log')
# ax1.grid()
# ax2.plot(P_u, thickness)
# ax2.set_title("Laminar flame thickness vs pressure.")
# ax2.set_ylabel("Laminar Flame Thickness [mm]")
# ax2.set_xlabel("Pressure [atm]")
# ax2.set_xscale('log')
# ax2.grid()
# fig.tight_layout()

plt.show()
