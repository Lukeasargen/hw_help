import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

def PM_function(M, gamma=1.4):
    part1 = np.sqrt( (gamma+1)/(gamma-1) )
    part2 = np.arctan( np.sqrt( (M*M - 1)*(gamma-1)/(gamma+1) ) )
    part3 = np.arctan( np.sqrt( M*M - 1) )
    return np.degrees(part1*part2 - part3)

def PM_function_from_nu(M_guess, nu, gamma=1.4):
    def f(M_guess, nu, gamma):
        return nu - PM_function(M_guess, gamma)
    M_guess = np.array([M_guess]*len(nu))
    M = fsolve(f, M_guess, args=(nu, gamma))
    return M

gamma = 1.4
M1 = 2.032
M2 = 2.460

# mu is mach wave angle relative to the flow
alpha1 = mu1 = np.degrees(np.arcsin(1/M1))
print(f"Leading Wave: {mu1 = } deg")
mu2 = np.degrees(np.arcsin(1/M2))
print(f"Trailing Wave: {mu2 = } deg")

nu1 = PM_function(M1, gamma=gamma)
print(f"{nu1 = } deg")
nu2 = PM_function(M2, gamma=gamma)
print(f"{nu2 = } deg")
delta = theta2 = nu2 - nu1
print(f"{theta2 = } deg")

alpha = mu1 - mu2 + theta2
print(f"Total fin angle: {alpha = } deg")

alpha2 = alpha1 - alpha
print(f"{alpha2 = } deg")
