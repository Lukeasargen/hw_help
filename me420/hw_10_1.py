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

def M_from_F(F, gamma=1.4):
    def f(M, F, gamma):
        return F - Fanno_from_M(M, gamma)
    M_guess = np.array([0.1]*len(F))
    M = fsolve(f, M_guess, args=(F, gamma))
    return M

# air constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

L = 50  # m, length
D = 0.065  # m, diameter
T1 = 240  # K, inlet temperature
P1 = 360*1000  # Pa, inlet pressure
M1 = 0.15  # inlet mach number
f = 0.0266  # Darcy friction factor

F = f*L/D
print(f"{F = }")
F1_star = Fanno_from_M(M1, gamma=gamma)
print(f"{F1_star = }")
F2_star = F1_star - F
print(f"{F2_star = }")
M2 = M_from_F([F2_star], gamma=gamma)[0]
print(f"{M2 = }")
L_star = F1_star*(D/f)
print(f"{L_star = } m")
