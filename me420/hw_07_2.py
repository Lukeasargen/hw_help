import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def shock_M2_from_M1(M1, gamma=1.4):
    top = 1 + ((gamma-1)/2)*M1**2
    bot = gamma*M1**2 - ((gamma-1)/2)
    return np.sqrt(top/bot)

def shock_P_ratio_from_M1(M1, gamma=1.4):  # P2/P1
    return(2*gamma*M1**2 - gamma + 1)/(gamma+1)

def shock_T_ratio_from_M1(M1, gamma=1.4):  # T2/T1
    part1 = 1 + 2*gamma*(M1*M1-1)/(gamma+1)
    part2 = (2 + (gamma-1)*M1*M1)/((gamma+1)*M1*M1)
    return part1*part2

# Constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

# Given: ambient
T_atm = 308.15  # K, ambient temperature
P_atm = 100.4e3  # Pa, ambient pressure
rho_atm = P_atm/(R*T_atm)  # kg/m^3, ambient density
print(f"{rho_atm = } kg/m^3")

# first atomic bomb measurements
E_1 = (18.6*1000)*(4184*1000)  # J, 18.6 kilotons of TNT to J
R_1 = 210  # m, radius at T_1
T_1 = 0.1  # s, time
print(f"{E_1 = :e} {E_1}")

# dimensionless_parameter
a = 1/5  # density
b = -1/5  # energy
c = -2/5  # time

pi_1 = R_1 * rho_atm**a * E_1**b * T_1**c
print(f"{pi_1 = }")

# bomb 2
R_2 = 482  # m, radius
T_2 = 0.113  # s, time
E_2 = (pi_1 * R_2**-1 * rho_atm**-a * T_2**-c)**(1/b)
print(f"{E_2 = :e} {E_2}")

# Plots are for bomb 2
t = np.linspace(0.1, 3.0, 100)
df = pd.DataFrame(t, columns=["t [s]"])
# Use similarity variable to find Radius with time
df["R [m]"] = pi_1 * rho_atm**-a * E_2**-b * df["t [s]"]**-c
K = pi_1 * rho_atm**-a * E_2**-b * -c
print(f"{K = }")
df["Vs [m/s]"] = K * df["t [s]"]**(-c-1)
a1 = np.sqrt(gamma*R*T_atm)
print(f"{a1 = } m/s")
df["Ms"] = df["Vs [m/s]"]/a1
# Static properties are the same for moving normal shocks
df["P2/P1"] = shock_P_ratio_from_M1(df["Ms"], gamma=gamma)
df["Pas [kPa]"] = df["P2/P1"]*P_atm/1000
df["T2/T1"] = shock_T_ratio_from_M1(df["Ms"], gamma=gamma)
df["Tas [K]"] = df["T2/T1"]*T_atm

# 2e) Find Ms=1 -> Vs=a1
t_mach = ((gamma*R*T_atm)/(K*K))**(1/(-2*c-2))
print(f"{t_mach = } s")
R_mach = pi_1 * rho_atm**-a * E_2**-b * t_mach**-c
print(f"{R_mach = } m")

df.loc[df["Ms"]<1, "Ms"] = 1
df["Mas"] = shock_M2_from_M1(df["Ms"], gamma=gamma)  # If Ms<1, set Ms=1

# 2f Find Pas=300 kPa
Pf = 300e3   # Pa, bc P_atm is in Pa
Pf_ratio = Pf/P_atm
print(f"{Pf_ratio = }")
Mf = np.sqrt((Pf_ratio*(gamma+1) + (gamma-1))/(2*gamma))
print(f"{Mf = }")
tf = ((Mf*a1)/K)**(1/(-c-1))
print(f"{tf = }")
Rf = pi_1 * rho_atm**-a * E_2**-b * tf**-c
print(f"{Rf = } m")

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Plot 2d
fig, (ax_R, ax_T, ax_P) = plt.subplots(nrows=3, ncols=1, figsize=(5,4))
plots = [
    [ax_R, "R [m]", "Radius"],
    [ax_T, "Tas [K]", "Temperature after shock"],
    [ax_P, "Pas [kPa]", "Pressure after shock"]
]
x = "t [s]"
for ax, y, y_label in plots:
    df.plot(x=x, y=y, ax=ax, legend=False)
    ax.set_title(f"{y_label} over time.")
    ax.set_ylabel(y)
    ax.set_xlabel("")
    ax.grid()
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontname('Arial')
        label.set_fontsize(SMALL_SIZE)
ax_P.set_xlabel(x)
fig.tight_layout()

# Plot 2e
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
df.plot(x=x, y="Ms", ax=ax, linestyle="-")
df.plot(x=x, y="Mas", ax=ax, linestyle=":")
ax.set_title("Mach number over time.")
ax.set_ylabel("M")
ax.grid()
fig.tight_layout()
plt.show()

