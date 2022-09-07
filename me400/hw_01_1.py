from re import M
import numpy as np
import matplotlib.pyplot as plt

# CO2
a = 22.26
b = 5.981e-2
c = -3.501e-5
d = 7.469e-9

# cp = a + bT + cT^2 + dT^3, KJ/Kmol*K
def f(T):
    return a + b*T + c*T**2 + d*T**3

# int cp = a T + 1/2 b T^2 + 1/3 c T^3 + 1/4 d T^4
def g(T):
    return a*T + (1/2)*b*T**2 + (1/3)*c*T**3 + (1/4)*d*T**4

M_CO2= 44.01 # g/mol, CO2 molar mass
T1 = 1000  # K
T2 = 353  # K
m_dot = 50  # kg/s

T = np.linspace(273, 1800, 100)
cp = f(T)/M_CO2
print(cp)

cp_int = g(T2)-g(T1)  # kJ/Kmol
print("Integral", cp_int)
W_cv = m_dot/M_CO2 * cp_int  # kWatts
print("W_cv", W_cv)

# Check part A
t_avg = 0.5*(T1+T2)  # K
print(f"{t_avg=}")
cp_avg = f(t_avg)/M_CO2  # KJ/Kg*K
print(f"{cp_avg=}")
W_cv = m_dot*cp_avg*(T1-T2)  # kWatts
print(f"{W_cv=}")
