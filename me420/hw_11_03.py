import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve


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

def beta_from_M_theta(beta_guess, M, theta, gamma=1.4):
    def f(beta_guess, M, theta, gamma):
        b_rad = np.radians(beta_guess)
        top = (M*np.sin(b_rad))**2 - 1
        bot = M*M*(gamma + np.cos(2*b_rad)) + 2
        return np.tan(np.radians(theta)) - 2*(top/bot)/np.tan(b_rad)
    M = np.array(M)
    beta_guess = np.array([beta_guess]*len(M))
    beta = fsolve(f, beta_guess, args=(M, theta, gamma))
    return beta

gamma = 1.4

M1 = 2.58  # mach number, inlet
T1 = 331.3  # K, inlet temperature
P1 = 192.2*1000  # Pa, inlet pressure
theta1 = 13.9  # degrees, wall angle

def oblique_shock_reflection(M1, T1, P1, theta, gamma=1.4):
    beta1 = beta_from_M_theta(beta_guess=10, M=[M1], theta=theta, gamma=gamma)[0]
    M1n = M1*np.sin(np.radians(beta1))
    M2n = shock_M2_from_M1(M1n, gamma=gamma)
    M2 = M2n/np.sin(np.radians(beta1-theta1))
    P2_P1_ratio = shock_P_ratio_from_M1(M1n, gamma=gamma)
    P2 = P2_P1_ratio*P1
    T2_T1_ratio = shock_T_ratio_from_M1(M1n, gamma=gamma)
    T2 = T2_T1_ratio*T1
    return M2, T2, P2, beta1


M2, T2, P2, beta1 = oblique_shock_reflection(M1, T1, P1, theta1, gamma=1.4)
print(f"{beta1 = } deg")
print(f"{M2 = }")
print(f"{P2/1000 = } kPa")
print(f"{T2 = } K")

theta2 = theta1
M3, T3, P3, beta2 = oblique_shock_reflection(M2, T2, P2, theta2, gamma=1.4)
print(f"{beta2 = } deg")
print(f"{M3 = }")
print(f"{P3/1000 = } kPa")
print(f"{T3 = } K")

phi = beta2-theta2
print(f"{phi = } deg")
