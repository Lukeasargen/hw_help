import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

def M_to_area_ratio(M, gamma=1.4):
    inv_mach_number = 1/M
    base_first_part = 2/(gamma+1)
    base_second_part = 1 + ((gamma-1)/2)*M**2
    exponent_t_ratio = (gamma+1)/(2*(gamma-1))
    return inv_mach_number * (base_first_part*base_second_part)**exponent_t_ratio

P0 = 975  # kPa
Patm = 20.3  # kPa
gamma = 1.36


M_squared = ( (P0/Patm)**((gamma-1)/gamma)-1 ) / ( (gamma-1)/2 )
print(f"{M_squared=}")
Me = np.sqrt(M_squared)
print(f"{Me=}")

area_ratio = M_to_area_ratio(Me, gamma)
print(f"{area_ratio=}")

