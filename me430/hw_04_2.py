import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

Ru = 8314.5  # J/kmol*K
P0 = 101325  # Pa

A = 2e12
b = 0
Ea = 2e8  # J/kmol
T = 2100  # K
Kp = 9.0

Ta = Ea/Ru  # K
print(f"{Ta=}")
Kf = A*(T**b)*np.exp(-Ta/T)
print(f"{Kf=}")
Kc = Kp * (P0)/(Ru*T)
print(f"{Kc=}")
Kr = Kf/Kc
print(f"{Kr=}")


