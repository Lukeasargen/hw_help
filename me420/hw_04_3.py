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

def M_to_P_ratio(M, gamma=1.4):
    # returns P/P0
    base = 1 + ((gamma-1)/2)*M**2
    exponent = -gamma/(gamma-1)
    return base**exponent

# Given
# Air
gamma = 1.4
R = 0.287  # kJ/Kg*K

P0 = 500  # kPa
T0 = 400  # K

x_r = np.array([
    [0.00, 0.3500],
    [0.02, 0.3209],
    [0.04, 0.2934],
    [0.06, 0.2673],
    [0.08, 0.2429],
    [0.10, 0.2202],
    [0.12, 0.1993],
    [0.14, 0.1803],
    [0.16, 0.1632],
    [0.18, 0.1483],
    [0.20, 0.1356],
    [0.22, 0.1254],
    [0.24, 0.1179],
    [0.26, 0.1131],
    [0.28, 0.1115]
])
df = pd.DataFrame(x_r, columns=["x [m]", "Radius [m]"])
df["Area [m^2]"] = np.pi*df["Radius [m]"]**2

# 3a) choked pressure = critical pressure
P_star = M_to_P_ratio(M=1, gamma=1.4)*P0
print(f"P* = {P_star}")

# 3b) Given Pb. Find Me
def P_to_M(P, P0, gamma=1.4):
    first_part = 2/(gamma-1)
    base = P0/P
    exponent = (gamma-1)/gamma
    return np.sqrt(first_part*((base**exponent) - 1))

# For finding the mach number from the area ratio
def f(M, area_ratio, gamma):
    return area_ratio - M_to_area_ratio(M, gamma)

# 3b, 3c, 3d, 4e
exit_idx = df["x [m]"].idxmax()  # exit is the highest x
A_exit = df.loc[exit_idx, "Area [m^2]"]
print(f"Ae = {A_exit} m^2")
print(f"Re = {df.loc[exit_idx, 'Radius [m]']} m")

fig, (ax_m, ax_p) = plt.subplots(nrows=1, ncols=2, figsize=(9,4.5))
x = "x [m]"
ax_m.set_title("Mach number as function of distance.",fontsize=10)
ax_m.set_xlabel(x)
ax_m.set_ylabel("Mach Number")
ax_p.set_title("Pressure as function of distance.",fontsize=10)
ax_p.set_xlabel(x)
ax_p.set_ylabel("Pressure [kPa]")

for Pb_label in [450, 400, 350, 300, 250, P_star]:
    Pb = max(P_star, Pb_label)
    M_exit = P_to_M(P=Pb, P0=P0, gamma=1.4)
    print(f"Me (Pb={Pb_label:.1f}) = {M_exit}")
    A_star = A_exit/M_to_area_ratio(M_exit, gamma=1.4)
    print(f"A* (Pb={Pb_label:.1f}) = {A_star} m^2")

    df[f"A/A* (Pb={Pb_label:.0f})"] = df["Area [m^2]"]/A_star
    x0_sub = np.array([0.01]*len(df))
    subsonic = fsolve(f, x0_sub, args=(df[f"A/A* (Pb={Pb_label:.0f})"], gamma))
    df[f"M (Pb={Pb_label:.0f})"] = subsonic
    df[f"P/P0 (Pb={Pb_label:.0f})"] = M_to_P_ratio(df[f"M (Pb={Pb_label:.0f})"], gamma)
    df[f"P (Pb={Pb_label:.0f})[kPa]"] = df[f"P/P0 (Pb={Pb_label:.0f})"]*P0

    label = f"{Pb_label:.0f} kPa"
    ax_m.plot(df[x], df[f"M (Pb={Pb_label:.0f})"], label=label)
    ax_m.annotate(label, (df[x].iloc[-1], df[f"M (Pb={Pb_label:.0f})"].iloc[-1]))
    ax_p.plot(df[x], df[f"P (Pb={Pb_label:.0f})[kPa]"], label=label)
    ax_p.annotate(label, (df[x].iloc[-1], df[f"P (Pb={Pb_label:.0f})[kPa]"].iloc[-1]))

ax_m.grid()
# ax_m.legend()
ax_p.grid()
# ax_p.legend()
fig.tight_layout()
plt.show()

