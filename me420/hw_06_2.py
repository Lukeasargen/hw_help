from cProfile import label
from traceback import print_tb
from matplotlib import markers
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

def shock_temp_ratio(M1, M2, gamma=1.4):
    top = 1 + ((gamma-1)/2)*M1**2
    bot = 1 + ((gamma-1)/2)*M2**2
    return top/bot

def fanno_curve(M1, M2, gamma=1.4):
    first_part = M1/M2
    return first_part*np.sqrt(shock_temp_ratio(M1, M2, gamma))

def rayleigh_curve(M1, M2, gamma=1.4):
    return (1+gamma*M1**2)/(1+gamma*M2**2)

def hugoniot_curve(density_ratio, gamma=1.4):
    gamma_ratio = (gamma+1)/(gamma-1)
    top = gamma_ratio*density_ratio - 1
    bot = gamma_ratio - density_ratio
    return top/bot

def M2_normal_shock_from_M1(M1, gamma=1.4):
    top = 1 + ((gamma-1)/2)*M1**2
    bot = gamma*M1**2 - ((gamma-1)/2)
    return np.sqrt(top/bot)

# Given
# Air
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1004.5  # J/Kg*K

M1 = 3.4
T1 = 258  # K
s1 = 1590  # J/Kg*K
h1 = cp*T1  # J/Kg
print(f"{h1 = } J/Kg")

# Find the properties after the shock
M2 = M2_normal_shock_from_M1(M1, gamma=gamma)
print(f"{M2 = }")
P_ratio  = fanno_curve(M1, M2, gamma=gamma)
print(f"{P_ratio = }")
T_ratio = P_ratio**2 * (M2/M1)**2
print(f"{T_ratio = }")
T2 = T_ratio*T1
print(f"{T2 = } K")
h2 = cp*T2
print(f"{h2 = } J/Kg")
s2 = cp*np.log(T_ratio) - R*np.log(P_ratio) + s1
print(f"{s2 = } J/Kg*K")

# PLOTS Corrected from hw_05
samples = 250  # Number of test points for M2

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))

# M2_list = np.linspace(0.1, 4.0, samples)
M2_list = np.exp(np.linspace(np.log(0.1), np.log(8.0), samples))
df = pd.DataFrame(M2_list, columns=["M2"])
df["M2/M1"] = df["M2"]/M1

# Steps
# 1) Find P2/P1. Use fanno, rayleigh, hugoniot
# 2) Ideal Gas + Speed of Sound -> T2/T1 = (P2/P1)^2 * (M2/M1)^2
# 3) Find T2 = T2/T1 * T1
# 4) Find h2 = cp*T2
# 5) Find s2 = cp*ln(T2/T1) - R*ln(P2/P1) + s1

# Step 1
df["Fanno P2/P1"] = fanno_curve(M1, df["M2"], gamma=gamma)
df["Rayleigh P2/P1"] = rayleigh_curve(M1, df["M2"], gamma=gamma)

# Create the density ratios
df["rho2/rho1"] = np.sqrt(1/shock_temp_ratio(M1, df["M2"], gamma=gamma))*(M1/df["M2"])
df["Hugoniot P2/P1"] = hugoniot_curve(df["rho2/rho1"], gamma=gamma)

# Steps 2-5
curve_names = ["Fanno", "Rayleigh", "Hugoniot"]
linestyles = ["-", "-.", "--"]
for curve_name, linestyle in zip(curve_names, linestyles):
    df[f"{curve_name} T2/T1"] = df[f"{curve_name} P2/P1"]**2 * df["M2/M1"]**2
    df[f"{curve_name} T2"] = df[f"{curve_name} T2/T1"]*T1
    df[f"{curve_name} h2"] = cp*df[f"{curve_name} T2"]
    df[f"{curve_name} s2"] = cp*np.log(df[f"{curve_name} T2/T1"]) - R*np.log(df[f"{curve_name} P2/P1"]) + s1
    df.plot(x=f"{curve_name} s2", y=f"{curve_name} h2", ax=ax, label=f"{curve_name}", linestyle=linestyle)

print(df)

ax.set_xlim([1300, 2500])
ax.set_ylim([1e5, 1e6])

ax.plot(s1, h1, marker="*")
ax.annotate(f"{M1=:.2f}", (s1, h1))
ax.plot(s2, h2, marker="*")
ax.annotate(f"{M2=:.2f}", (s2, h2))
ax.set_title("Mollier diagram")
ax.set_xlabel("Entropy [J/Kg*K]")
ax.set_ylabel("Enthalpy [J/Kg]")
fig.tight_layout()
plt.show()
