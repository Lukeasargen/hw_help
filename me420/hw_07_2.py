import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

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

def shock_a_ratio_from_M1(M1, gamma=1.4):  # a2/a1
    top = 2*(gamma-1)*(gamma*M1*M1 - 1/(M1*M1) - (gamma-1))
    bot = (gamma+1)**2
    return np.sqrt(1 + top/bot)

def shock_rho_ratio_from_M1(M1, gamma=1.4):  # rho2/rho1
    term = 2*(1 - 1/(M1*M1)) / (gamma+1)
    return 1/(1-term)

def Mas_from_M1(M1, gamma=1.4):
    return M1/np.sqrt(shock_T_ratio_from_M1(M1, gamma=gamma)) - shock_M2_from_M1(M1, gamma=gamma)

def M1_from_Mas(Mas, gamma=1.4):
    def f(M, Mas, gamma):
        return Mas - Mas_from_M1(M, gamma)
    M_guess = np.array([1.0]*len(Mas))
    M1 = fsolve(f, M_guess, args=(Mas, gamma))
    return M1

# Constants
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K
J_per_kiloton = 4.184e12  # Joules in kiloton of TNT
J_per_megaton = 4.184e15  # Joules in megaton of TNT

# Given: ambient
T_atm = 308.15  # K, ambient temperature
P_atm = 100.4e3  # Pa, ambient pressure
rho_atm = P_atm/(R*T_atm)  # kg/m^3, ambient density
print(f"{rho_atm = } kg/m^3")

# first atomic bomb measurements
kiloton_1 = 18.6
R_1 = 210  # m, radius at T_1
T_1 = 0.1  # s, time
E_1 = kiloton_1*J_per_kiloton  # J, 18.6 kilotons of TNT to J
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
megaton_2 = E_2/J_per_megaton
print(f"{megaton_2 = :e} {megaton_2}")

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

# Define the shock to be Ms=1 after Ms hits 1
df.loc[df["Ms"]<1, "Ms"] = 1  # If Ms<1, set Ms=1

# Static properties are the same for normal shocks
df["P2/P1"] = shock_P_ratio_from_M1(df["Ms"], gamma=gamma)
df["Pas [kPa]"] = df["P2/P1"]*P_atm/1000
df["T2/T1"] = shock_T_ratio_from_M1(df["Ms"], gamma=gamma)
df["Tas [K]"] = df["T2/T1"]*T_atm
df["a2"] = np.sqrt(gamma*R*df["Tas [K]"])
# df["a2"] = np.sqrt(df["T2/T1"])*a1

# Use M2 to find Mas
df["M2"] = shock_M2_from_M1(df["Ms"], gamma=gamma)
# df["V2 [m/s]"] = df["M2"]*df["a2"]
# df["Vas [m/s]"] = df["Vs [m/s]"]-df["V2 [m/s]"]
# df["Mas"] = df["Vas [m/s]"]/df["a2"]
# df["Mas"] = df["Vs [m/s]"]/df["a2"] - df["M2"]
df["Mas"] = Mas_from_M1(df["Ms"], gamma=gamma)
df.loc[df["Mas"]<0, "Mas"] = 0  # If Ms<1, set Ms=1

# 2e) Find Ms=1 -> Vs=a1
t_e1 = ((gamma*R*T_atm)/(K*K))**(1/(-2*c-2))
print(f"{t_e1 = } s")
R_e1 = pi_1 * rho_atm**-a * E_2**-b * t_e1**-c
print(f"{R_e1 = } m")

# 2e) Find Mas=1
M1_e2 = M1_from_Mas([1], gamma=gamma)[0]
print(f"{M1_e2 = }")
t_e2 = ((M1_e2*a1)/K)**(1/(-c-1))
print(f"{t_e2 = } s")
R_e2 = pi_1 * rho_atm**-a * E_2**-b * t_e2**-c
print(f"{R_e2 = } m")

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

