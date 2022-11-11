import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve


def T0_T0star_ratio_from_M(M, gamma=1.4):  # T01/T0*
    top = (2+(gamma-1)*M*M)*(1+gamma)*M*M
    bot = (1+gamma*M*M)**2
    return top/bot

def M_from_T0_T0star_ratio(ratio, gamma=1.4):
    def f(M, ratio, gamma):
        return ratio - T0_T0star_ratio_from_M(M, gamma)
    M_guess = np.array([0.5]*len(ratio))
    M = fsolve(f, M_guess, args=(ratio, gamma))
    return M

def T_ratio_from_M1_M2(M1, M2, gamma=1.4):  # T2/T1
    top = M2*(1+gamma*M1*M1)
    bot = M1*(1+gamma*M2*M2)
    return (top/bot)**2

def P_ratio_from_M1_M2(M1, M2, gamma=1.4):  # P2/P1
    top = 1+gamma*M1*M1
    bot = 1+gamma*M2*M2
    return top/bot

# air constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

diameter = 0.135  # m, diameter
T1 = 505  # K, inlet temperature
P1 = 360*1000  # Pa, inlet pressure
V1 = 110  # m/s, inlet velocity
Q = 5406*1000  # W, heat added

area = 0.25*np.pi*diameter**2
print(f"{area = } m^2")

rho1 = P1/(R*T1)
print(f"{rho1 = } kg/m^3")
a1 = np.sqrt(gamma*R*T1)
print(f"{a1 = } m/s")
M1 = V1/a1
print(f"{M1 = }")
T01 = T1 + (V1**2)/(2*cp)
print(f"{T01 = } K")

m_dot = rho1*area*V1
print(f"{m_dot = } kg/s")
q = Q/m_dot 
print(f"{q = } kJ/kg")
T02 = T01 + q/cp
print(f"{T02 = } K")

T01_T0star_ratio = T0_T0star_ratio_from_M(M1, gamma=gamma)
print(f"{T01_T0star_ratio = }")
T02_T0star_ratio = (T02/T01)*T01_T0star_ratio
# print(f"{T02_T0star_ratio = }")
M2 = M_from_T0_T0star_ratio([T02_T0star_ratio], gamma=gamma)[0]
print(f"{M2 = }")

T0star = T01/T01_T0star_ratio
print(f"{T0star = } K")
Q_max = m_dot*cp*(T0star-T01)
print(f"{Q_max/1000 = } kW")
print("Not choked." if Q<Q_max else "CHOKED.")

T2_T1_ratio = T_ratio_from_M1_M2(M1, M2, gamma=gamma)
# print(f"{T2_T1_ratio = }")
T2 = T2_T1_ratio*T1
print(f"{T2 = } K")
P2_P1_ratio = P_ratio_from_M1_M2(M1, M2, gamma=gamma)
# print(f"{P2_P1_ratio = }")
P2 = P2_P1_ratio*P1
print(f"{P2/1000 = } kPa")
rho2 = P2/(R*T2)
print(f"{rho2 = } kg/m^3")
a2 = np.sqrt(gamma*R*T2)
print(f"{a2 = } m/s")
V2 = a2*M2
print(f"{V2 = } m/s")
