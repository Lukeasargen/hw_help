from re import sub
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def M_to_area_ratio(M, gamma=1.4):
    inv_mach_number = 1/M
    base_first_part = 2/(gamma+1)
    base_second_part = 1 + ((gamma-1)/2)*M**2
    exponent_t_ratio = (gamma+1)/(2*(gamma-1))
    return inv_mach_number * (base_first_part*base_second_part)**exponent_t_ratio

# Given
gamma = 1.4
P0 = 220  # kPa
T0 = 300  # K

# Calculate area along the nozzle from class
# x = np.linspace(0.0, 0.3, 13)
# A = np.array([
#     0.031415927,
#     0.031151243,
#     0.0303351,
#     0.029093379,
#     0.02755741,
#     0.025855012,
#     0.02410492,
#     0.022413343,
#     0.02087258,
#     0.019561706,
#     0.01854931,
#     0.01789831,
#     0.017671459
# ])
# df = pd.DataFrame(x, columns=["x [m]"])
# df["Area [m^2]"] = A
# exit_row_idx = df["x [m]"].idxmax()
# A_exit = df.loc[exit_row_idx, "Area [m^2]"]
# print(f"A_exit: {A_exit} m^2")
# df["A/A*"] = df["Area [m^2]"]/A_exit

# homework problem
df = pd.DataFrame([0.0], columns=["x [m]"])
df["Area [m^2]"] = 0.015
A_exit = 0.00445
df["A/A*"] = df["Area [m^2]"]/A_exit


# a) False Position Method
# Goal is to iterate the mach number to match the area ratio
def false_position_method(guess_bounds, area_ratio, gamma=1.4, max_steps=100):
    assert len(guess_bounds)==2
    # Make initial guesses, (rows, 2), min and max guess
    mach_est = np.tile(np.array(guess_bounds), (len(area_ratio),1))
    for i in range(max_steps):
        # Uses guess to find area ratio
        ar = M_to_area_ratio(M=mach_est, gamma=gamma)
        # Linear interpolation
        new_M = mach_est[:,0]+(mach_est[:,1]-mach_est[:,0])*(area_ratio-ar[:,0])/(ar[:,1]-ar[:,0])
        # Smooth the step size for better stability
        alpha = 0.4
        new_M = alpha*new_M + (1-alpha)*mach_est[:,0]
        # Update the guesses
        mach_est[:,0] = mach_est[:,1]
        mach_est[:,1] = new_M
        # If the guesses are very close, stop iterations
        max_diff = np.max(np.abs(mach_est[:,0]-mach_est[:,1]))
        if max_diff<1e-6: break
    # print(f"Steps {i+1}")
    return mach_est[:,1]  # Return the latest guess
        
subsonic = false_position_method([0.1, 0.9], area_ratio=df["A/A*"], gamma=gamma)
supersonic = false_position_method([2.0, 3.0], area_ratio=df["A/A*"], gamma=gamma)
df["M FP sub"] = subsonic
df["M FP sup"] = supersonic

# b) Newton's method
# Goal is to iterate the gradient to 0
def newton_method(guess, area_ratio, gamma=1.4, max_steps=100):
    mach_est = np.array([guess]*len(area_ratio))
    a = 2/(gamma+1)
    b = (gamma-1)/2
    c = (gamma+1)/(2*(gamma-1))
    for i in range(max_steps):
        G = area_ratio - M_to_area_ratio(M=mach_est, gamma=gamma)
        base = (a + a*b*mach_est**2)
        first_part = (-1/(mach_est**2)) * (base**c)
        second_part = (1/mach_est)*(c*2*a*b*mach_est) * base**(c-1)
        dG = - first_part - second_part
        alpha = 1.1  # scale the step size a little when g/dG is very small
        mach_est = mach_est - alpha*G/dG
        if np.abs(G/dG).max() < 1e-6: break
    # print(f"Steps {i+1}")
    return mach_est

subsonic = newton_method(0.2, area_ratio=df["A/A*"], gamma=gamma)
supersonic = newton_method(1.5, area_ratio=df["A/A*"], gamma=gamma)
df["M Newton sub"] = subsonic
df["M Newton sup"] = supersonic

# 3) Use existing program
from scipy.optimize import fsolve
def f(M, area_ratio, gamma):
    return area_ratio - M_to_area_ratio(M, gamma)
x0_sub = np.array([0.2]*len(df))
x0_sup = np.array([1.5]*len(df))
subsonic = fsolve(f, x0_sub, args=(df["A/A*"], gamma))
supersonic = fsolve(f, x0_sup, args=(df["A/A*"], gamma))
df["M fsolve sub"] = subsonic
df["M fsolve sup"] = supersonic

print(df)
