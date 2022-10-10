from cProfile import label
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

def rayleigh_curve(M1, M2, gamma):
    return (1+gamma*M1**2)/(1+gamma*M2**2)

def M2_from_M1(M1, gamma=1.4):
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

h1 = cp*T1

# 2d)
M2 = M2_from_M1(M1, gamma=1.4)
print(f"{M2 = }")
T_ratio = shock_temp_ratio(M1, M2, gamma=gamma)
print(f"{T_ratio = }")
P_ratio  = fanno_curve(M1, M2, gamma=gamma)
print(f"{P_ratio = }")
T2 = T_ratio*T1
print(f"{T2 = }")
h2 = cp*T2
print(f"{h2 = }")
s2 = cp*np.log(T_ratio) - R*np.log(P_ratio) + s1
print(f"{s2 = }")


# PLOT

M2_list = np.exp(np.linspace(np.log(0.1), np.log(5.0), 50))
print()
df = pd.DataFrame(M2_list, columns=["M2"])

df["T2/T1"] = shock_temp_ratio(M1, df["M2"], gamma=gamma)
df["h2"] = cp*df["T2/T1"]*T1

df["Fanno P2/P1"] = fanno_curve(M1, df["M2"], gamma=gamma)
fanno_delta_entropy = cp*np.log(df["T2/T1"]) - R*np.log(df["Fanno P2/P1"])
df["Fanno s2"] = fanno_delta_entropy + s1

df["Rayleigh P2/P1"] = rayleigh_curve(M1, df["M2"], gamma=gamma)
rayleigh_delta_entropy = cp*np.log(df["T2/T1"]) - R*np.log(df["Rayleigh P2/P1"])
df["Rayleigh s2"] = rayleigh_delta_entropy + s1

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))

df.plot(x="Fanno s2", y="h2", ax=ax, label="Fanno")
df.plot(x="Rayleigh s2", y="h2", ax=ax, label="Rayleigh")

ax.plot(s1, h1, marker="*")
ax.annotate(f"{M1=:.2f}", (s1, h1))

ax.plot(s2, h2, marker="*")
ax.annotate(f"{M2=:.2f}", (s2, h2))

ax.set_title("Mollier diagram")
ax.set_xlabel("Entropy [J/Kg*K]")
ax.set_ylabel("Enthalpy [J/Kg]")

fig.tight_layout()
plt.show()
