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

def rayleigh_curve(M1, M2, gamma=1.4):
    return (1+gamma*M1**2)/(1+gamma*M2**2)

# air constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

diameter = 0.135  # m, diameter
T1 = 505  # K, inlet temperature
P1 = 360*1000  # Pa, inlet pressure
V1 = 110  # m/s, inlet velocity

samples = 100
s1 = 2850  # J/kg*K

# Does the rayleigh from the first part
Q = np.linspace(0, 6132, samples) * 1000  # W, heat added
area = 0.25*np.pi*diameter**2
rho1 = P1/(R*T1)
a1 = np.sqrt(gamma*R*T1)
M1 = V1/a1
T01 = T1 + (V1**2)/(2*cp)
m_dot = rho1*area*V1
q = Q/m_dot 
T02 = T01 + q/cp
T01_T0star_ratio = T0_T0star_ratio_from_M(M1, gamma=gamma)
T02_T0star_ratio = (T02/T01)*T01_T0star_ratio
M2 = M_from_T0_T0star_ratio(T02_T0star_ratio, gamma=gamma)
T0star = T01/T01_T0star_ratio
Q_max = m_dot*cp*(T0star-T01)
T2_T1_ratio = T_ratio_from_M1_M2(M1, M2, gamma=gamma)
T2 = T2_T1_ratio*T1
P2_P1_ratio = P_ratio_from_M1_M2(M1, M2, gamma=gamma)
P2 = P2_P1_ratio*P1
rho2 = P2/(R*T2)
a2 = np.sqrt(gamma*R*T2)
V2 = a2*M2
s2 = cp*np.log(T2_T1_ratio) - R*np.log(P2_P1_ratio) + s1

# theoretical rayleigh line
M2_rayleigh = np.linspace(M1, 1.0, samples)
P1_P2_rayleigh = rayleigh_curve(M1, M2_rayleigh, gamma=1.4)
T2_T1_rayleigh = P1_P2_rayleigh**2 * (M2_rayleigh/M1)**2
T2_rayleigh = T2_T1_rayleigh*T1
s2_rayleigh = cp*np.log(T2_T1_rayleigh) - R*np.log(P1_P2_rayleigh) + s1

# Plot this line
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
ax.plot(s2, T2, label="Rayleigh part 1", marker="o")
ax.plot(s2_rayleigh, T2_rayleigh, label="Rayleigh theoretical")
ax.set_title("T-s diagram")
ax.set_xlabel("Entropy [J/Kg*K]")
ax.set_ylabel("Temperature [K]")
ax.grid()
ax.legend()
fig.tight_layout()
plt.show()
