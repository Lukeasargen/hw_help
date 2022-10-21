import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# Constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

# altitude of 37k ft
T1 = 216.650  # K, temperature
P1 = 21.6627  # kPa, pressure
rho1 = P1/(R*T1)  # kg/m^3, ambient density
print(f"{rho1 = } kg/m^3")

P0 = 52.3  # kPa, pitot pressire, stagnation pressure

# Incompressible flow bernoulli
M1 = np.sqrt((2*(P0-P1))/(gamma*R*T1*rho1))
print(f"{M1 = } incompressible bernoulli")
# Compressible flow, isentropic relation
M1 = (2/(gamma-1))*((P0/P1)**((gamma-1)/gamma)-1)
print(f"{M1 = } compressible subsonic")

def M_to_P1_P02_ratio(M, gamma=1.4):
    base1 = (1 - gamma + 2*gamma*M*M)/(gamma+1)
    part1 = (base1)**(1/(gamma-1))
    base2 = ((gamma+1)/2)*M*M
    part2 = base2**(-gamma/(gamma-1))
    return part1*part2

def f(M, P1_P02_ratio, gamma):
    return P1_P02_ratio - M_to_P1_P02_ratio(M, gamma)

P1_P02_ratio = P1/P0
print(f"{P1_P02_ratio = }")
M_guess = np.array([1.5])
M1 = fsolve(f, M_guess, args=(P1_P02_ratio, gamma))[0]
print(f"{M1 = } compressible supersonic")

a1 = np.sqrt(gamma*R*T1)
print(f"{a1 = } m/s")

V1 = a1*M1
print(f"{V1 = } m/s")
