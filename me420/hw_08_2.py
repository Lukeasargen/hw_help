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

def Vcs_from_Ms(Ms, gamma1, a1):
    return (2/(gamma1+1))*(Ms-(1/Ms))*a1

def Ms_from_Vcs(Vcs, gamma1, a1):
    def f(M, Vcs, gamma1, a1):
        return Vcs - Vcs_from_Ms(M, gamma1, a1)
    M_guess = np.array([1.5]*len(Vcs))
    Ms = fsolve(f, M_guess, args=(Vcs, gamma1, a1))
    return Ms

# Low pressure side, Argon
P1 = 101.3e3  # Pa
T1 = 320.15  # K
gamma1 = 1.667
R1 = 208.0  # J/Kg*K

# High pressure side, Air
P4 = 1250e3  # Pa
T4 = 500  # K
gamma4 = 1.4
R4 = 287.0  # J/Kg*K

# assume diaphragm is at x=0
xL = -12  # m, low pressure side
xR = 18  # m, high pressure side

# Calculations
def shock_tube(low_pressure_state, high_pressure_state, x):
    P1, T1, gamma1, R1 = low_pressure_state
    P4, T4, gamma4, R4 = high_pressure_state
    xL, xR = x
    a1 = np.sqrt(gamma1*R1*T1)
    a4 = np.sqrt(gamma4*R4*T4)
    rho1 = P1/(R1*T1)
    rho4 = P4/(R4*T4)
    P4_P1_ratio = P4/P1
    Ms = Ms_from_P4_P1([P4_P1_ratio], gamma=(gamma1, gamma4), a=(a1, a4))[0]
    Vs = Ms*a1
    P2_P1_ratio = shock_P_ratio_from_M1(Ms, gamma=gamma1)
    P2 = P2_P1_ratio*P1
    T2_T1_ratio = shock_T_ratio_from_M1(Ms, gamma=gamma1)
    T2 = T2_T1_ratio*T1
    a2 = np.sqrt(gamma1*R1*T2)
    V3 = Vcs = Vas = V2 = Vcs_from_Ms(Ms, gamma1, a1)
    M2 = V2/a2
    a3_a4_ratio = 1 - ((gamma4-1)/2)*(V3/a4)
    a3 = a3_a4_ratio*a4
    T3_T4_ratio = (a3/a4)**2
    T3 = T3_T4_ratio*T4
    P3_P4_ratio = (a3/a4)**((2*gamma4)/(gamma4-1))
    P3 = P3_P4_ratio*P4
    rho3_rho4_ratio = (a3/a4)**((2)/(gamma4-1))
    rho3 = rho3_rho4_ratio*rho4
    tR = abs(xR/Vs)
    tL = abs(xL/a4)
    # Reflection xR
    Msr = Msr_from_Ms(np.array([Ms]), a1, a2)[0]
    Vsr = Vsr_from_Ms(Ms, gamma1, a1)
    M5 = shock_M2_from_M1(Msr, gamma=gamma1)
    T5_T2_ratio = shock_T_ratio_from_M1(Msr, gamma=gamma1)
    T5 = T5_T2_ratio*T2
    P5_P2_ratio = shock_P_ratio_from_M1(Msr, gamma=gamma1)
    P5 = P5_P2_ratio*P2

    states =  {
        "P1": P1, "T1": T1, "gamma1": gamma1, "R1": R1,
        "P4": P4, "T4": T4, "gamma4": gamma4, "R4": R4,
        "xL": xL, "xR": xR,
        "a1": a1, "a4": a4,
        "rho1": rho1, "rho4": rho4,
        "P4_P1_ratio": P4_P1_ratio,
        "Ms": Ms, "Vs": Vs,
        "P2_P1_ratio": P2_P1_ratio, "P2": P2,
        "T2_T1_ratio": T2_T1_ratio, "T2": T2,
        "V3": V3, "Vcs": Vcs, "Vas": Vas, "V2": V2,
        "a2": a2, "M2": M2,
        "a3_a4_ratio": a3_a4_ratio, "a3": a3,
        "T3_T4_ratio": T3_T4_ratio, "T3": T3,
        "P3_P4_ratio": P3_P4_ratio, "P3": P3,
        "rho3_rho4_ratio": rho3_rho4_ratio, "rho3": rho3,
        "Vcs-a3": Vcs-a3,
        "tR": tR, "tL": tL,
        "Msr": Msr, "Vsr": Vsr,
        "M5": M5,
        "T5_T2_ratio": T5_T2_ratio, "T5": T5,
        "P5_P2_ratio": P5_P2_ratio, "P5": P5,
        "gamma4": gamma4,
    }
    return states

states1 = shock_tube(
    low_pressure_state = (P1, T1, gamma1, R1),
    high_pressure_state = (P4, T4, gamma4, R4),
    x = (xL, xR),
)
print(f"{states1=}")

# 2d)
V3d = states1["a4"]/(1+((gamma4-1)/2))
print(f"{V3d = } m/s")
Msd = Ms_from_Vcs([V3d], gamma1, states1["a1"])[0]
print(f"{Msd = }")
P4_P1_ratio_d = P4_P1_from_Ms(Msd, gamma=(gamma1, gamma4), a=(states1["a1"], states1["a4"]))
P4d = P4_P1_ratio_d*P1
print(f"{P4_P1_ratio_d = }")
print(f"{P4d/1000 = } kPa")

states2 = shock_tube(
    low_pressure_state = (P1, T1, gamma1, R1),
    high_pressure_state = (P4d, T4, gamma4, R4),
    x = (xL, xR),
)
print(f"{states2=}")

# Plot states1 on x-t diagram
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,3))
ax.plot([0,xR],[0,states1["tR"]], label="Vs", linestyle="-")
ax.plot([0,xL],[0,states1["tL"]], label="a4", linestyle=":")
ax.plot([0,states1["tR"]*states1["Vcs"]],[0,states1["tR"]], label="Vcs", linestyle="--")
ax.plot([0,states1["tR"]*states1["Vcs-a3"]],[0,states1["tR"]], label="Vcs-a3", linestyle="-.")

ax.set_xlabel("x [m]")
ax.set_ylabel("t [s]")
ax.set_ylim(bottom=0)
ax.set_xlim([xL, xR])
ax.grid()
ax.legend()
fig.tight_layout()
plt.show()
