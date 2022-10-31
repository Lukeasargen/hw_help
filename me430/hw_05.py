import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint

# Constants
C_kg_kmol = 12
H_kg_kmol = 1
O_kg_kmol = 16
N_kg_kmol = 14
Ru = 8315  # J/Kmol*K

# Inputs
# Fuel, 1 molar unit, CxHyOz
x = 1; y = 4; z = 0
equiv_ratio = 0.8  # <1.0
P = 5*101325  # atm -> Pa, reactants
T = 1000  # K, reactants

# Calculations below
inv_equiv_ratio = 1/equiv_ratio
O2_vol = 0.21
N2_vol = 0.79  # includes the other intert gases like Argon
N2_mole_per_O2 = 3.76  # N2_vol/O2_vol 3.76 moles N2 per 1 mole O2
# molecular weights
fuel_kg_kmol = x*C_kg_kmol + y*H_kg_kmol + z*O_kg_kmol
O2_kg_kmol = 2*O_kg_kmol
N2_kg_kmol = 2*N_kg_kmol
CO2_kg_kmol = C_kg_kmol + 2*O_kg_kmol
H2O_kg_kmol = 2*H_kg_kmol + O_kg_kmol

# Global one-step mechanism
print("CxHyOz + a*(O2+3.76N2) -> b*CO2 + c*H2O + e*O2 + f*N2")
a = x + y/4 - z/2
b = x
c = y/2
e = a
f = N2_mole_per_O2*a
print(f"{a=}")
print(f"{inv_equiv_ratio*a=}")
print(f"{b=}")
print(f"{c=}")
print(f"{e=}")
print(f"{(inv_equiv_ratio-1)*e=}")
print(f"{f=}")
print(f"{inv_equiv_ratio*f=}")

reactants = {
    "Gas": ["CxHyOz", "O2", "N2"],
    "kmol": [1, inv_equiv_ratio*a*1, inv_equiv_ratio*a*N2_mole_per_O2],
    "kg/kmol": [fuel_kg_kmol, O2_kg_kmol, N2_kg_kmol],
}
reactants = pd.DataFrame(reactants)

products = {
    "Gas": ["CO2", "H2O", "O2", "N2"],
    "kmol": [b, c, (inv_equiv_ratio-1)*e, inv_equiv_ratio*f],
    "kg/kmol": [CO2_kg_kmol, H2O_kg_kmol, O2_kg_kmol, N2_kg_kmol],
}
products = pd.DataFrame(products)

for df, name in [(reactants, "Reactants"), (products, "Products")]:
    total_kmols = df["kmol"].sum()
    print(f"{name} Total kmol: {total_kmols} kmol")
    df["Mole Fraction"] = df["kmol"]/total_kmols
    total_molecular_weight = (df["Mole Fraction"]*df["kg/kmol"]).sum()
    print(f"{name} Total Molecular Weight: {total_molecular_weight} kg/kmol")
    print(f"{name} Total Mass: {total_kmols*total_molecular_weight} kg")
    df["Mass Fraction"] = (df["Mole Fraction"]*df["kg/kmol"])/total_molecular_weight

print("INTIAL REACTANTS.")
reactants["Molar Concentration"] = (reactants["Mole Fraction"]*P)/(Ru*T)
print(reactants)

print("GLOBAL ONE-STEP APPROXIMATIONS PRODCUTS.")
print("CONSTANT-TEMPERATURE and CONSTANT-PRESSURE PRODUCTS")
products["Molar Concentration"] = (products["Mole Fraction"]*P)/(Ru*T)
print(products)


# Find the intitial-time-rate-of-change of fuel molar concentration
print("initial time-rate-of-change of fuel molar concentration")
# Values are given, they also come from Table 5.1
A = 1.3e8  # (kmol/m3)^(1-m-n) / s 
b = 0
m = -0.3  # order methane
n = 1.3  # order oxygen
EA = 202525*1000  # J/Kmol, given in kJ but converted to J

def get_value(df, gas, col):
    return float(df[df["Gas"]==gas][col])

fuel_c = get_value(reactants, "CxHyOz", "Molar Concentration")
O2_c = get_value(reactants, "O2", "Molar Concentration")
fuel_k = -A * T**b * np.exp(-EA/(Ru*T)) * fuel_c**m * O2_c**n
print(f"d[CxHyOz]/dt = {fuel_k} kmol/m3-s")  # negative because it's depleted
X_fuel = get_value(reactants, "CxHyOz", "Mole Fraction")
X_O2 = get_value(reactants, "O2", "Mole Fraction")
fuel_k = -A * T**b *np.exp(-EA/(Ru*T)) * X_fuel**m * X_O2**n * (P/(Ru*T))**(m+n)
print(f"d[CxHyOz]/dt = {fuel_k} kmol/m3-s")  # negative because it's depleted

def ODE_rates(molar_concentration, t, T, A, b, m, n, EA):  # t is required for ODE solver
    fuel_c, O2_c, CO2_c, H2O_c = molar_concentration
    d_fuel = -A * T**b * np.exp(-EA/(Ru*T)) * fuel_c**m * O2_c**n
    d_O2 = 2*d_fuel  # there are 2.5 kmol reac - 0.5 kmol prod = 2 kmol net O2
    d_C02 = -d_fuel
    d_H2O = -2*d_fuel
    omega_dot = np.array([d_fuel, d_O2, d_C02, d_H2O])
    sum_c = np.sum(molar_concentration)
    sum_omega = np.sum(omega_dot)
    return omega_dot # - molar_concentration*(sum_omega/sum_c)

# constant-pressure, constant-temperature
# dP/dt = dT/dt = 0
molar_concentration = np.array([
    get_value(reactants, "CxHyOz", "Molar Concentration"),
    get_value(reactants, "O2", "Molar Concentration"),
    0,  # CO2
    0,  # H2O
])  # kmol/m^3, Initial conditions
print(f"{molar_concentration = }")

t = np.linspace(0, 150, 10000)
sol = odeint(ODE_rates, molar_concentration, t, args=(T, A, b, m, n, EA))

idxs = (sol[:,0]>molar_concentration[0]*0.001)  # CH4 0.1% of initial
sol = sol[idxs]
t = t[idxs]

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))
ax.plot(t, sol[:, 0], label='CH4', linestyle="solid", color='red')
ax.plot(t, sol[:, 1], label='O2', linestyle="dashed", color='blue')
ax.plot(t, sol[:, 2], label='CO2', linestyle="dashdot", color='green')
ax.plot(t, sol[:, 3], label='H2O', linestyle="dotted", color='black')
ax.legend(loc='best')
ax.set_title("Time evolution of molar concentrations")
ax.set_xlabel('Time [s]')
ax.set_ylabel("Molar Concentration [kmol/m^3]")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.grid()
fig.tight_layout()
plt.show()
