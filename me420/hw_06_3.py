import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

def M_to_area_ratio(M, gamma=1.4):
    inv_mach_number = 1/M
    base_first_part = 2/(gamma+1)
    base_second_part = 1 + ((gamma-1)/2)*M**2
    exponent_t_ratio = (gamma+1)/(2*(gamma-1))
    return inv_mach_number * (base_first_part*base_second_part)**exponent_t_ratio

def M_to_P_ratio(M, gamma=1.4):  # P/P0
    base = 1 + ((gamma-1)/2)*M**2
    exponent = -gamma/(gamma-1)
    return base**exponent

def P_to_M(P, P0, gamma=1.4):
    first_part = 2/(gamma-1)
    base = P0/P
    exponent = (gamma-1)/gamma
    return np.sqrt(first_part*((base**exponent) - 1))

def M2_from_M1(M1, gamma=1.4):
    top = 1 + ((gamma-1)/2)*M1**2
    bot = gamma*M1**2 - ((gamma-1)/2)
    return np.sqrt(top/bot)

def shock_P_ratio_from_M1(M1, gamma=1.4):  # P2/P1
    return(2*gamma*M1**2 - gamma + 1)/(gamma+1)

def shock_T_ratio_from_M(M1, M2, gamma=1.4):  # T2/T1
    top = 1 + ((gamma-1)/2)*M1**2
    bot = 1 + ((gamma-1)/2)*M2**2
    return top/bot

def shock_P0_ratio_from_M(M1, M2, gamma=1.4):  # P02/P01
    base = 1/shock_T_ratio_from_M(M1, M2, gamma)
    exponent = (gamma+1)/(2*(gamma-1))
    return (M1/M2) * base**exponent

def shock_P02_ratio_from_M1(M1, gamma=1.4):  # P1/P02
    first_part = shock_P_ratio_from_M1(M1, gamma=gamma)**(1/(gamma-1))
    second_part = (((gamma+1)/2) * M1**2)**(-gamma/(gamma-1))
    return first_part*second_part

def shock_rho_ratio_from_M1(M1, gamma=1.4):  # rho2/rho1
    return ((gamma+1)*M1**2) / (2 + (gamma-1)*M1**2)

# For finding the mach number from the area ratio
def f(M, area_ratio, gamma):
    return area_ratio - M_to_area_ratio(M, gamma)

# Given
# Air
gamma = 1.4
R = 0.287  # kJ/Kg*K
cp = 1005  # J/Kg*K

P01 = 765  # kPa
T01 = 483  # K

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

# Geometry at the throat
throat_idx = df["Area [m^2]"].idxmin()  # throat is the smallest area
A_throat = df.loc[throat_idx, "Area [m^2]"]
x_throat = df.loc[throat_idx, "x [m]"]
print(f"{x_throat = } m")
print(f"{A_throat = } m^2")
# You can also calculate P1*
P1_star = M_to_P_ratio(M=1, gamma=gamma)*P01
print(f"P1* = {P1_star}")
# Find back pressure for subsonic and supersonic divergence with no shocks
exit_idx = df["x [m]"].idxmax()  # exit is the highest x
A_exit = df.loc[exit_idx, "Area [m^2]"]
print(f"{A_exit = } m^2")
area_ratio_exit = A_exit/A_throat  # assume the throat is choked
print(f"{area_ratio_exit = }")
M_subsonic = fsolve(f, 0.5, args=(area_ratio_exit, gamma))[0]
print(f"{M_subsonic = }")
P_exit_subsonic = M_to_P_ratio(M=M_subsonic, gamma=gamma)*P01
print(f"{P_exit_subsonic = }")
M_supersonic = fsolve(f, 1.5, args=(area_ratio_exit, gamma))[0]
print(f"{M_supersonic = }")
P_exit_supersonic = M_to_P_ratio(M=M_supersonic, gamma=gamma)*P01
print(f"{P_exit_supersonic = }")

print(f"Shock occurs if the back pressure is between {P_exit_subsonic:.2f} kPa and {P_exit_supersonic:.2f} kPa")

# Case where the shock position is known
shock_position_x = 0.9

# The throat must be choked
A1_star = A_throat
print(f"{A1_star = }")
df["A/A1*"] = df["Area [m^2]"]/A1_star

# Calculate the mach number for converging nozzle
x0_sub = np.array([0.1]*len(df))
subsonic = fsolve(f, x0_sub, args=(df["A/A1*"], gamma))
df["M"] = subsonic

# Replace the mach number for the supersonic divergence
x0_sup = np.array([1.5]*len(df))
supersonic = fsolve(f, x0_sup, args=(df["A/A1*"], gamma))
df.loc[throat_idx:, "M"] = supersonic[throat_idx:]
df.loc[throat_idx, "M"] = 1

# Calculate the properties if there was no shock
temp_ratio = 1 / (1 + ((gamma-1)/2)*df[f"M"]**2)  # T/T0
df["T"] = temp_ratio * T01
pressure_exponent = gamma/(gamma-1)
df["P"] = P01 * temp_ratio**pressure_exponent  # P01 * P/P01 = P01 * (T/T01)^(gamma/gamma-1)

# Find the area ratio at the shock position and M1 before the shock
shock_area_ratio = interp1d(df["x [m]"], df["A/A1*"])(shock_position_x)
shock_M1 = interp1d(df["x [m]"], df["M"])(shock_position_x)
shock_P1 = interp1d(df["x [m]"], df["P"])(shock_position_x)
print(f"{shock_area_ratio = }")
print(f"{shock_M1 = }")
print(f"{shock_P1 = }")

# Find the region after the shock using the position index
# This is the index right after the shock, because the shock position is assumed to have supersonic conditions
shock_idx = df[df['x [m]'].gt(shock_position_x)].index[0]

# Calculate properties accross the shock
shock_M2 = M2_from_M1(shock_M1, gamma=1.4)
print(f"{shock_M2 = }")
shock_T_ratio = shock_T_ratio_from_M(shock_M1, shock_M2, gamma=gamma)  # T2/T1
print(f"{shock_T_ratio = }")
shock_P_ratio = shock_P_ratio_from_M1(shock_M1, gamma=gamma)  # P2/P1
print(f"{shock_P_ratio = }")
shock_P0_ratio = shock_P0_ratio_from_M(shock_M1, shock_M2, gamma=gamma)  # P02/P01
print(f"{shock_P0_ratio = }")
shock_P02_ratio = shock_P02_ratio_from_M1(shock_M1, gamma=gamma)  # P1/P02
print(f"{shock_P02_ratio = }")
shock_rho_ratio = shock_rho_ratio_from_M1(shock_M1, gamma=gamma)  # rho2/rho1
print(f"{shock_rho_ratio = }")
shock_P2 = shock_P1*shock_P_ratio
print(f"{shock_P2 = }")

# Calculate new critical conditions
P02 = shock_P0_ratio*P01
print(f"{P02 = }")
P2_star = M_to_P_ratio(M=1, gamma=gamma)*P02
print(f"P2* = {P2_star}")
A2_star = A1_star/shock_P0_ratio
print(f"A2* = {A2_star}")
df["A/A2*"] = df["Area [m^2]"]/A2_star

# Replace the mach number with subsonic case after the shock
x0_sub = np.array([0.1]*len(df))
subsonic = fsolve(f, x0_sub, args=(df[f"A/A2*"], gamma))
df.loc[shock_idx:, "M"] = subsonic[shock_idx:]

# Re-calculate the properties after the shock
temp_ratio = 1 / (1 + ((gamma-1)/2)*df[f"M"]**2)  # T/T0
pressure_exponent = gamma/(gamma-1)
df.loc[shock_idx:, "P"] = P02 * temp_ratio[shock_idx:]**pressure_exponent  # P01 * P/P01 = P01 * (T/T01)^(gamma/gamma-1)

# Back pressure at exit
back_pressure = df.iloc[-1]["P"]
print(f"{back_pressure = }")

fig, (ax_m, ax_p) = plt.subplots(nrows=2, ncols=1, figsize=(5,6))
x = "x [m]"
ax_m.set_title("Mach number as function of distance.")
ax_m.set_xlabel(x)
ax_m.set_ylabel("Mach Number")
df.plot(x=x, y="M", ax=ax_m)

ax_p.set_title("Pressure as function of distance.")
ax_p.set_xlabel(x)
ax_p.set_ylabel("Pressure [kPa]")
df.plot(x=x, y="P", ax=ax_p)

ax_m.grid()
ax_p.grid()
fig.tight_layout()
plt.show()
