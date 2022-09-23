import numpy as np
P0 = 165  # kPa
T0 = 493.5  # K

# a) Find T* and P* for air
gamma = 1.4  # air
T_star_air = T0 * 2 / (gamma+1)
print(f"{T_star_air=}")
P_star_air = P0 * (2/(gamma+1))**((gamma-1)/gamma)
print(f"{P_star_air=}")

# b) Find T* and P* for helium
gamma = 1.66
T_star_H = T0 * 2 / (gamma+1)
print(f"{T_star_H=}")
P_star_H = P0 * (2/(gamma+1))**((gamma-1)/gamma)
print(f"{P_star_H=}")

# c) Both are choked because Pb < P_star
Pb = 99.65  # kPa

# d) masss flow rate
Ae = 0.013  # m^2
Me = 1

def mass_flow_rate(P0, T0, A, M, gamma=1.4, Ri=287.0):
    """ P0=Pa, T0=K, A=m^2, M=mach number,
        gamma=heat ratio, Ri=species gas constant J/kg*K"""
    constant = np.sqrt(gamma/Ri)
    stagnation = P0*A*np.sqrt(1/T0)
    base = 1 + ((gamma-1)/2)*M**2
    exponent = -((gamma+1)/(2*(gamma-1)))
    mach_effects = M*(base**exponent)
    return mach_effects * constant * stagnation

# Air
m_dot_air = mass_flow_rate(P0=P0*1e3, T0=T0, A=Ae, M=Me,
                        gamma=1.4, Ri=287.0)
print(f"{m_dot_air=}")

m_dot_H = mass_flow_rate(P0=P0*1e3, T0=T0, A=Ae, M=Me,
                        gamma=1.66, Ri=2077.1)
print(f"{m_dot_H=}")
