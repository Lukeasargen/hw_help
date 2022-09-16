import matplotlib.pyplot as plt
import numpy as np

def mass_flow_rate(P0, T0, A, M, gamma=1.4, Ri=287.0):
    """ P0=Pa, T0=K, A=m^2, M=mach number,
        gamma=heat ratio, Ri=species gas constant J/kg*K"""
    constant = np.sqrt(gamma/Ri)
    stagnation = P0*A*np.sqrt(1/T0)
    base = 1 + ((gamma-1)/2)*M**2
    exponent = -((gamma+1)/(2*(gamma-1)))
    mach_effects = M*(base**exponent)
    return mach_effects * constant * stagnation

kg_s = mass_flow_rate(
            P0=600000,
            T0=500,
            A=0.0079485,
            M=1.0,
            gamma=1.4,
            Ri=287.0
            )

print(kg_s)

gamma = 1.4
M = 1
Ri = 287.0
constant = np.sqrt(gamma/Ri)
print(f"{constant=}")
base = 1 + ((gamma-1)/2)*M**2
exponent = -((gamma+1)/(2*(gamma-1)))
mach_effects = M*(base**exponent)
print(f"{mach_effects=}")

print(f"{1000*constant*mach_effects=}")

