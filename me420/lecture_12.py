import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def M_to_area_ratio(M, gamma=1.4):
    inv_mach_number = 1/M
    base_first_part = 2/(gamma+1)
    base_second_part = 1 + ((gamma-1)/2)*M**2
    exponent_t_ratio = (gamma+1)/(2*(gamma-1))
    return inv_mach_number * (base_first_part*base_second_part)**exponent_t_ratio

def M_to_P_ratio(M, gamma=1.4):
    base = 1 + ((gamma-1)/2)*M**2
    exponent = -gamma/(gamma-1)
    return base**exponent

# Given
gamma = 1.4
P0 = 220  # kPa
T0 = 300  # K
Pb = 50.0  # kPa

x_a = np.array([
    [0, 0.031415927],
    [0.025, 0.031151243],
    [0.05, 0.0303351],
    [0.075, 0.029093379],
    [0.1, 0.02755741],
    [0.125, 0.025855012],
    [0.15, 0.02410492],
    [0.175, 0.022413343],
    [0.2, 0.02087258],
    [0.225, 0.019561706],
    [0.25, 0.01854931],
    [0.275, 0.01789831],
    [0.3, 0.017671459],
    [0.325, 0.01807656],
    [0.35, 0.01925422],
    [0.375, 0.021124432],
    [0.4, 0.02360668],
    [0.425, 0.026600575],
    [0.45, 0.02997014],
    [0.475, 0.033534303],
    [0.5, 0.03706572],
    [0.525, 0.040297166],
    [0.55, 0.04293412],
    [0.575, 0.044674717],
    [0.6, 0.045238934]
])
df = pd.DataFrame(x_a, columns=["x [m]", "Area [m^2]"])
df["Radius [m]"] = np.sqrt(df["Area [m^2]"]/np.pi)

throat_idx = df["Area [m^2]"].idxmin()  # throat is the smallest area
A_throat = df.loc[throat_idx, "Area [m^2]"]
print(f"A_throat = {A_throat} m^2")
df["A/A*"] = df["Area [m^2]"]/A_throat

exit_idx = df["x [m]"].idxmax()  # exit is the highest x
A_exit = df.loc[exit_idx, "Area [m^2]"]
print(f"A_exit = {A_exit} m^2")

P_star = Pb/P0
print(f"P* = {P_star}")

from scipy.optimize import fsolve
def f(M, area_ratio, gamma):
    return area_ratio - M_to_area_ratio(M, gamma)
x0_sub = np.array([0.2]*len(df))
x0_sup = np.array([1.5]*len(df))
subsonic = fsolve(f, x0_sub, args=(df["A/A*"], gamma))
supersonic = fsolve(f, x0_sup, args=(df["A/A*"], gamma))
df["M"] = subsonic
df.loc[throat_idx:, "M"] = supersonic[throat_idx:]
df.loc[throat_idx, "M"] = 1

df["P/P0"] = M_to_P_ratio(df["M"], gamma)
df["P [kPa]"] = df["P/P0"]*P0

print(df)

x = "x [m]"
y_left = "M"
y_right = "P/P0"

fig = plt.figure()
ax_left = fig.add_subplot(1, 1, 1)
ax_left.set_title("CD Nozzle",fontsize=10)
ax_left.set_xlabel("x [m]")

ax_left.plot(df[x], df[y_left])
ax_left.set_ylabel(y_left)
ax_left.set_ylim(bottom=0)

ax_right = ax_left.twinx()
ax_right.plot(df[x], df[y_right])
ax_right.set_ylabel(y_right)
ax_right.set_ybound([0, 1])

ax_left.set_xbound([min(df[x]), max(df[x])])
ax_left.grid()
fig.tight_layout()
# plt.show()
