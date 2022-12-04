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

L = 154.1  # m, length
D = 0.075  # m, diameter
T1 = 350  # K, inlet temperature
P1 = 240*1000  # Pa, inlet pressure
V1 = 55.5  # inlet velocity
eta = 0.0008e-3  # m, average roughness
dynamic_viscosity1 = 18.03e-6  # Ns/m^2, at 350K

def fanno_pipe(L, D, T1, P1, V1, eta, dynamic_viscosity):
    states = {
        "L": L, "Dh": D,
        "P1": P1, "T1": T1, "V1": V1, "dynamic_viscosity": dynamic_viscosity,
    }
    states["rho1"] = P1/(R*T1)
    states["a1"] = np.sqrt(gamma*R*T1)
    states["M1"] = states["V1"]/states["a1"]
    states["Re1"] = (states["rho1"]*states["V1"]*states["Dh"])/states["dynamic_viscosity"]
    states["f1"] = churchill_equation(states["Re1"], states["Dh"], eta)
    states["F"] = states["f1"]*states["L"]/states["Dh"]
    states["F1_star"] = Fanno_from_M(states["M1"], gamma=gamma)
    states["F2_star"] = states["F1_star"] - states["F"]
    states["M2"] = M_from_F(states["F2_star"], gamma=gamma)
    states["L_star"] = states["F1_star"]*(states["Dh"]/states["f1"])
    states["T2_T1_ratio"] = T_ratio_from_M1_M2(states["M1"], states["M2"], gamma=gamma)
    states["T2"] = states["T2_T1_ratio"]*states["T1"]
    states["a2"] = np.sqrt(gamma*R*states["T2"])
    states["V2"] = states["a2"]*states["M2"]
    states["rho2"] = states["rho1"]*states["V1"]/states["V2"]
    states["P2"] = states["rho2"]*R*states["T2"]
    states["P2_P1_ratio"] = states["P2"]/states["P1"]
    states["ds"] = cp*np.log(states["T2_T1_ratio"]) - R*np.log(states["P2_P1_ratio"])
    return states

L = np.array([10, 50, 100, L])
states1 = fanno_pipe(L, D, T1, P1, V1, eta, dynamic_viscosity1)
for k,v in states1.items():
    print(k, v)

# c)
dynamic_viscosity2 = 20.27e-6  # Ns/m^2, at 339K
Re2 = (states1["rho2"]*states1["V2"]*states1["Dh"])/dynamic_viscosity2
print(f"{Re2 = }")
f2 = churchill_equation(Re2, states1["Dh"], eta)
print(f"{f2 = }")

# d)
samples = 50
L = np.linspace(0, states1["L_star"], samples)
states = fanno_pipe(L, D, T1, P1, V1, eta, dynamic_viscosity1)

# e)
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
# ax.plot(states["L"], states["M2"])
# ax.set_xlabel("L [m]")
# ax.set_ylabel("M2 []")
# ax.grid()
# fig.tight_layout()
# plt.show()

s1 = 2500
s2 = s1 + states["ds"]

def shock_temp_ratio(M1, M2, gamma=1.4):
    top = 1 + ((gamma-1)/2)*M1**2
    bot = 1 + ((gamma-1)/2)*M2**2
    return top/bot

def fanno_curve(M1, M2, gamma=1.4):
    first_part = M1/M2
    return first_part*np.sqrt(shock_temp_ratio(M1, M2, gamma))

fanno_P2_P1_ratio = fanno_curve(states["M1"], states["M2"], gamma=gamma)
fanno_T2_T1_ratio = fanno_P2_P1_ratio**2 * (states["M2"]/states["M1"])**2
fanno_T2 = fanno_T2_T1_ratio*states["T1"]
fanno_s2 = cp*np.log(fanno_T2_T1_ratio) - R*np.log(fanno_P2_P1_ratio) + s1


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
ax.plot(states["T2"], s2, label="Part d", marker="o")
ax.plot(fanno_T2, fanno_s2, label="Fanno theoretical")
ax.set_xlabel("T2 [K]")
ax.set_ylabel("S [J/Kg*K]")
ax.grid()
ax.legend()
fig.tight_layout()

plt.show()
