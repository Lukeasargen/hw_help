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

def P_to_M(P, P0, gamma=1.4):
    first_part = 2/(gamma-1)
    base = P0/P
    exponent = (gamma-1)/gamma
    return np.sqrt(first_part*((base**exponent) - 1))

# For finding the mach number from the area ratio
def f(M, area_ratio, gamma):
    return area_ratio - M_to_area_ratio(M, gamma)

# Given
# Air
gamma = 1.4
R = 0.287  # kJ/Kg*K

P0 = 765  # kPa
T0 = 483  # K

x_r = np.array([
    [0.00, 0.3500],
    [0.05, 0.3303],
    [0.10, 0.3114],
    [0.15, 0.2936],
    [0.20, 0.2771],
    [0.25, 0.2622],
    [0.30, 0.2492],
    [0.35, 0.2382],
    [0.40, 0.2293],
    [0.45, 0.2229],
    [0.50, 0.2189],
    [0.53, 0.2179],
    [0.55, 0.2174],
    [0.58, 0.2177],
    [0.60, 0.2185],
    [0.65, 0.2222],
    [0.70, 0.2282],
    [0.75, 0.2366],
    [0.80, 0.2472],
    [0.85, 0.2596],
    [0.90, 0.2737],
    [0.95, 0.2890],
    [1.00, 0.3050],
    [1.05, 0.3213],
    [1.10, 0.3372],
    [1.15, 0.3521],
    [1.20, 0.3651],
    [1.25, 0.3755],
    [1.30, 0.3823]
])
df = pd.DataFrame(x_r, columns=["x [m]", "Radius [m]"])
df["Area [m^2]"] = np.pi*df["Radius [m]"]**2

throat_idx = df["Area [m^2]"].idxmin()  # throat is the smallest area
A_throat = df.loc[throat_idx, "Area [m^2]"]
x_throat = df.loc[throat_idx, "x [m]"]
print(f"{x_throat = } m")
print(f"{A_throat = } m^2")
exit_idx = df["x [m]"].idxmax()  # exit is the highest x
A_exit = df.loc[exit_idx, "Area [m^2]"]
print(f"{A_exit = } m^2")
P_star = M_to_P_ratio(M=1, gamma=gamma)*P0
print(f"P* = {P_star}")

# 2a) choked, subsonic diffuser
A_star = A_throat  # assume the throat is choked
area_ratio_exit = A_exit/A_star
print(f"{area_ratio_exit = }")

print("Part 2a")
x0_sub = np.array([0.5])
M_subsonic = fsolve(f, x0_sub, args=(area_ratio_exit, gamma))[0]
print(f"{M_subsonic = }")
P_ratio_exit_subsonic = M_to_P_ratio(M=M_subsonic, gamma=gamma)
print(f"{P_ratio_exit_subsonic = }")
P_exit_subsonic = P_ratio_exit_subsonic*P0
print(f"{P_exit_subsonic = }")

# 2c) choked, supersonic nozzle
print("Part 2c")
x0_sup = np.array([1.5])
M_supersonic = fsolve(f, x0_sup, args=(area_ratio_exit, gamma))[0]
print(f"{M_supersonic = }")
P_ratio_exit_supersonic = M_to_P_ratio(M=M_supersonic, gamma=gamma)
print(f"{P_ratio_exit_supersonic = }")
P_exit_supersonic = P_ratio_exit_supersonic*P0
print(f"{P_exit_supersonic = }")

# 2b, 2d
fontsize = 8
fig, ((ax_m, ax_p), (ax_T, ax_rho)) = plt.subplots(nrows=2, ncols=2, figsize=(8,4))
x = "x [m]"
ax_m.set_title("Mach number as function of distance.", fontsize=fontsize)
ax_m.set_xlabel(x, fontsize=fontsize)
ax_m.set_ylabel("Mach Number", fontsize=fontsize)
ax_p.set_title("Pressure as function of distance.", fontsize=fontsize)
ax_p.set_xlabel(x, fontsize=fontsize)
ax_p.set_ylabel("Pressure [kPa]", fontsize=fontsize)
ax_T.set_title("Temperature as function of distance.", fontsize=fontsize)
ax_T.set_xlabel(x, fontsize=fontsize)
ax_T.set_ylabel("Temperature [K]", fontsize=fontsize)
ax_rho.set_title("Density as function of distance.", fontsize=fontsize)
ax_rho.set_xlabel(x, fontsize=fontsize)
ax_rho.set_ylabel("Density [kg/m^3]", fontsize=fontsize)

conditions = [
    [P_exit_subsonic, "2a. choked subsonic"],
    [P_exit_supersonic, "2c. choked supersonic"],
]

for Pb, str_label in conditions:
    choked = Pb < P_exit_subsonic
    print(str_label, Pb)
    print(f"{choked=}")
    A_star = A_throat
    # Use back pressure to find star conditions ???
    # M_exit = P_to_M(P=Pb, P0=P0, gamma=gamma)
    # print(f"Me ({Pb=:.1f}) = {M_exit}")
    # area_ratio = M_to_area_ratio(M_exit, gamma=gamma)
    # print(f"aera_ratio ({Pb=:.1f}) = {area_ratio}")
    # A_star = A_throat/area_ratio
    # print(f"A* ({Pb=:.1f}) = {A_star} m^2")

    # Use star conditions to find mach number
    df[f"A/A* ({Pb=:.0f})"] = df["Area [m^2]"]/A_star
    x0_sub = np.array([0.1]*len(df))
    subsonic = fsolve(f, x0_sub, args=(df[f"A/A* ({Pb=:.0f})"], gamma))
    df[f"M ({Pb=:.0f})"] = subsonic
    if choked:  # Replace the solution for the diverging section
        x0_sup = np.array([1.5]*len(df))
        supersonic = fsolve(f, x0_sup, args=(df[f"A/A* ({Pb=:.0f})"], gamma))
        df.loc[throat_idx:, f"M ({Pb=:.0f})"] = supersonic[throat_idx:]
        df.loc[throat_idx, f"M ({Pb=:.0f})"] = 1

    temp_ratio = 1 + ((gamma-1)/2)*df[f"M ({Pb=:.0f})"]**2
    df[f"T ({Pb=:.0f})"] = T0/temp_ratio
    pressure_exponent = -gamma/(gamma-1)
    df[f"P ({Pb=:.0f})"] = P0 * temp_ratio**pressure_exponent
    density_exponent = -1/(gamma-1)
    df[f"p ({Pb=:.0f})"] = P0 * temp_ratio**density_exponent

    ax_m.plot(df[x], df[f"M ({Pb=:.0f})"], label=str_label)
    ax_p.plot(df[x], df[f"P ({Pb=:.0f})"], label=str_label)
    ax_T.plot(df[x], df[f"T ({Pb=:.0f})"], label=str_label)
    ax_rho.plot(df[x], df[f"p ({Pb=:.0f})"], label=str_label)

print(df)

ax_m.legend()
ax_m.grid()
ax_p.grid()
ax_T.grid()
ax_rho.grid()

fig.tight_layout()
plt.show()
