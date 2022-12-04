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

def churchill_equation(Re, Dh, eta):
    A = (-2.457*np.log((7/Re)**0.9+(0.27*eta/Dh)))**16
    B = (37530/Re)**16
    f = 8*((8/Re)**12 + (A+B)**-1.5)**(1/12)
    return f

def T_ratio_from_M1_M2(M1, M2, gamma=1.4):  # T2/T1
    top = 1 + ((gamma-1)/2)*M1**2
    bot = 1 + ((gamma-1)/2)*M2**2
    return top/bot


# air constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K
dynamic_viscosity = 24.99e-6  # Ns/m^2, at 450 K

L = 27.0  # m, length
D = 0.05  # m, diameter
T1 = 450  # K, inlet temperature
P1 = 220*1000  # Pa, inlet pressure
V1 = 85  # inlet velocity
eta = 0.08e-3  # m, average roughness

# a) churchill equation
Dh = D  # m, hydraulic diameter
rho1 = P1/(R*T1)
print(f"{rho1 = } kg/m^3")
a1 = np.sqrt(gamma*R*T1)
print(f"{a1 = } m/s")
M1 = V1/a1
print(f"{M1 = }")
Re = (rho1*V1*Dh)/dynamic_viscosity
print(f"{Re = }")
f1 = churchill_equation(Re, Dh, eta)
print(f"{f1 = }")

# b)
F = f1*L/D
print(f"{F = }")
F1_star = Fanno_from_M(M1, gamma=gamma)
print(f"{F1_star = }")
F2_star = F1_star - F
print(f"{F2_star = }")
M2 = M_from_F([F2_star], gamma=gamma)[0]
print(f"{M2 = }")
L_star = F1_star*(D/f1)

print(f"{L_star = } m")
T2_T1_ratio = T_ratio_from_M1_M2(M1, M2, gamma=gamma)
print(f"{T2_T1_ratio = }")
T2 = T2_T1_ratio*T1
print(f"{T2 = } K")
a2 = np.sqrt(gamma*R*T2)
print(f"{a2 = } m/s")
V2 = a2*M2
print(f"{V2 = } m/s")
rho2 = rho1*V1/V2
print(f"{rho2 = } kg/m^3")
P2 = rho2*R*T2
print(f"{P2/1000 = } kPa")

