import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

def shock_temp_ratio(M1, M2, gamma=1.4):  # T2/T1
    top = 1 + ((gamma-1)/2)*M1**2
    bot = 1 + ((gamma-1)/2)*M2**2
    return top/bot

def shock_P_ratio_from_M1(M1, gamma=1.4):  # P2/P1
    return(2*gamma*M1**2 - gamma + 1)/(gamma+1)

def shock_T_ratio_from_M1(M1, gamma=1.4):  # T2/T1
    part1 = 1 + 2*gamma*(M1*M1-1)/(gamma+1)
    part2 = (2 + (gamma-1)*M1*M1)/((gamma+1)*M1*M1)
    return part1*part2

def shock_M2_from_M1(M1, gamma=1.4):
    top = 1 + ((gamma-1)/2)*M1**2
    bot = gamma*M1**2 - ((gamma-1)/2)
    return np.sqrt(top/bot)

def P4_P1_from_Ms(Ms, gamma, a):
    g1, g4 = gamma
    a1, a4 = a
    exponent = (-2*g4)/(g4-1)
    base = 1 - ((g4-1)/(g1+1))*(a1/a4)*(Ms-(1/Ms))
    P_ratio = shock_P_ratio_from_M1(Ms, gamma=g1)
    return P_ratio*(base**exponent)

def Ms_from_P4_P1(P4_P1_ratio, gamma, a):
    def f(M, P4_P1_ratio, gamma, a):
        return P4_P1_ratio - P4_P1_from_Ms(M, gamma, a)
    M_guess = np.array([1.5]*len(P4_P1_ratio))
    Ms = fsolve(f, M_guess, args=(P4_P1_ratio, gamma, a))
    return Ms

def Msr_from_Ms(Ms, a1, a2):
    def f(M, Ms, a1, a2):
        return (M*M-1)/M - (a1/a2)*(Ms*Ms-1)/Ms
    M_guess = np.array([1.5]*len(Ms))
    Msr = fsolve(f, M_guess, args=(Ms, a1, a2))
    return Msr

def Vsr_from_Ms(Ms, gamma1, a1):
    k = (gamma1-1)/(gamma1+1)
    Vsr = a1*(1+2*k*(Ms*Ms-1))/Ms
    return Vsr

# Low pressure side, Air
P1 = 100e3  # Pa
T1 = 298.15  # K
gamma1 = 1.4
R1 = 287.0  # J/Kg*K

# High pressure side, Helium
P4 = 800e3  # Pa
T4 = 298.15  # K
gamma4 = 1.663
R4 = 2077.15  # J/Kg*K

# assume diaphragm is at x=0
xL = -15  # m, low pressure side
xR = 10  # m, high pressure side

# Calculations
a1 = np.sqrt(gamma1*R1*T1)
a4 = np.sqrt(gamma4*R4*T4)
print(f"{a1 = } m/s")
print(f"{a4 = } m/s")

rho1 = P1/(R1*T1)
rho4 = P4/(R4*T4)
print(f"{rho1 = } kg/m^2")
print(f"{rho4 = } kg/m^2")

P4_P1_ratio = P4/P1
print(f"{P4_P1_ratio = }")
Ms = M1 = Ms_from_P4_P1([P4_P1_ratio], gamma=(gamma1, gamma4), a=(a1, a4))[0]
Vs = Ms*a1
print(f"{Ms = }")
print(f"{Vs = } m/s")

# 2a) P2, P3, T2, T3, Ms, Vs, Vcv
P2_P1_ratio = shock_P_ratio_from_M1(Ms, gamma=gamma1)
P2 = P2_P1_ratio*P1
print(f"{P2_P1_ratio = }")
print(f"{P2/1000 = } kPa")

T2_T1_ratio = shock_T_ratio_from_M1(Ms, gamma=gamma1)
T2 = T2_T1_ratio*T1
print(f"{T2_T1_ratio = }")
print(f"{T2 = } K")

a2 = np.sqrt(gamma1*R1*T2)
V3 = Vcs = Vas = V2 = (2/(gamma1+1))*(Ms-(1/Ms))*a1
M2 = V2/a2
print(f"{a2 = } m/s")
print(f"{Vcs = } m/s")
print(f"{M2 = }")

a3_a4_ratio = 1 - ((gamma4-1)/2)*(V3/a4)
a3 = a3_a4_ratio*a4
print(f"{a3_a4_ratio = }")
print(f"{a3 = } m/s")

T3_T4_ratio = (a3/a4)**2
T3 = T3_T4_ratio*T4
print(f"{T3_T4_ratio = }")
print(f"{T3 = } K")

P3_P4_ratio = (a3/a4)**((2*gamma4)/(gamma4-1))
P3 = P3_P4_ratio*P4
print(f"{P3_P4_ratio = }")
print(f"{P3/1000 = } kPa")

rho3_rho4_ratio = (a3/a4)**((2)/(gamma4-1))
rho3 = rho3_rho4_ratio*rho4
print(f"{rho3_rho4_ratio = }")
print(f"{rho3 = } kg/m^2")

# 2b)
print(f"{Vs = } m/s")
print(f"{Vcs = } m/s")
print(f"{-a3+Vcs = } m/s")
print(f"{-a4 = } m/s")
tR = xR/Vs
print(f"{tR = } s")
tL = -xL/a4
print(f"{tL = } s")

# 2c) 
Msr = Msr_from_Ms(np.array([Ms]), a1, a2)[0]
print(f"{Msr = }")
Vsr = Vsr_from_Ms(Ms, gamma1, a1)
print(f"{Vsr = } m/s")
M5 = shock_M2_from_M1(Msr, gamma=gamma1)
print(f"{M5 = }")
T5_T2_ratio = shock_temp_ratio(Msr, M5, gamma=gamma1)
T5 = T5_T2_ratio*T2
print(f"{T5_T2_ratio = }")
print(f"{T5 = } K")
P5_P2_ratio = shock_P_ratio_from_M1(Msr, gamma=gamma1)
P5 = P5_P2_ratio*P2
print(f"{P5_P2_ratio = }")
print(f"{P5/1000 = } kPa")

# Plot these states on x-t diagram
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,3))
ax.plot([0,xR],[0,tR], label="Vs", linestyle="-")
ax.plot([0,xL],[0,tL], label="a4", linestyle=":")
ax.plot([0,tR*Vcs],[0,tR], label="Vcs", linestyle="--")
ax.plot([0,tR*(-a3+Vcs)],[0,tR], label="a3-Vcs", linestyle="-.")

ax.set_xlabel("x [m]")
ax.set_ylabel("t [s]")
ax.set_ylim(bottom=0)
ax.set_xlim([xL, xR])
ax.grid()
ax.legend()
fig.tight_layout()
plt.show()
