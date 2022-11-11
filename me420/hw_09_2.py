import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

def Fanno_from_M(M, gamma=1.4):
    part1 = (1-M*M)/(gamma*M*M)
    top = (gamma+1)*M*M
    bot = 2 + (gamma-1)*M*M
    part2 = (gamma+1)/(2*gamma)*np.log(top/bot)
    return part1 + part2

# air constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

L = 50  # m, length
D = 0.065  # m, diameter
T1 = 240  # K, inlet temperature
P1 = 360*1000  # Pa, inlet pressure
M1 = 0.13  # inlet mach number
f = 0.0216  # Darcy friction factor

F = f*L/D
print(f"{F = }")
F1_star = Fanno_from_M(M1, gamma=gamma)
print(f"{F1_star = }")
F2_star = F1_star - F
print(f"{F2_star = }")
L_star = F1_star*(D/f)
print(f"{L_star = } m")

samples = 100
M = np.linspace(0.2, 6, samples)
F = Fanno_from_M(M, gamma=gamma)
print(Fanno_from_M(M=100, gamma=gamma))
lim = ((gamma+1)*np.log((gamma+1)/(gamma-1)))/(2*gamma) - 1/gamma
print(lim)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
ax.plot(M, F)
ax.set_xlabel("M")
ax.set_ylabel("fL*/D")
ax.grid()
fig.tight_layout()
plt.show()
