import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def get_value(df, gas, col):
    return float(df[df["Gas"]==gas][col])

# Constants
C_kg_kmol = 12
H_kg_kmol = 1
O_kg_kmol = 16
N_kg_kmol = 14
Ru = 8315  # J/Kmol*K, reactants

# Inputs
# Fuel, 1 molar unit, CxHyOz
x = 3; y = 8; z = 0
equiv_ratio = 1.0  # <1.0
P_reac = 1*101325  # atm -> Pa, reactants
T_reac = 300  # K, reactants, unburned

# Table B.1, constant pressure adiabatic flame temperature for fuel-air
T_prod = 2267  # K
T_avg = (T_reac+T_prod)/2
print(f"{T_avg = } K")
# Average temperature in the heat release zone
T_avg_reaction = (T_avg+T_prod)/2
print(f"{T_avg_reaction = } K")

# Use the average temperautre to get cp_air and k_air
# Table C.1, properties of air at T_avg
cp_avg = 1189  # J/kg*K, specific heat
k_avg = 82e-3  # W/m*K, thermal conductivity

# Table 5.1, global one step reaction rates
# A = 8.6e11  # (kmol/m3)^(1-m-n) / s 
A = 4.836e9  # not in Table 5.1, but given in Example 8.2 for propane-air
b_coeff = 0
m = 0.1  # order fuel
n = 1.65  # order oxygen
TA = 15098  # K, activation temperature

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
print(f"{a = }")
print(f"{inv_equiv_ratio*a = }")
print(f"{b = }")
print(f"{c = }")
print(f"{e = }")
print(f"{(inv_equiv_ratio-1)*e = }")
print(f"{f = }")
print(f"{inv_equiv_ratio*f = }")

print("INTIAL REACTANTS.")
reactants = {
    "Gas": ["CxHyOz", "O2", "N2"],
    "kmol": [1, inv_equiv_ratio*a*1, inv_equiv_ratio*a*N2_mole_per_O2],
    "kg/kmol": [fuel_kg_kmol, O2_kg_kmol, N2_kg_kmol],
}
reactants = pd.DataFrame(reactants)
reactants_kmols = reactants["kmol"].sum()
print(f"Reactants Total kmol: {reactants_kmols} kmol")
reactants["Mole Fraction"] = reactants["kmol"]/reactants_kmols
reactants_molecular_weight = (reactants["Mole Fraction"]*reactants["kg/kmol"]).sum()
print(f"Reactants Total Molecular Weight: {reactants_molecular_weight} kg/kmol")
print(f"Reactants Total Mass: {reactants_kmols*reactants_molecular_weight} kg")
reactants["Mass Fraction"] = (reactants["Mole Fraction"]*reactants["kg/kmol"])/reactants_molecular_weight
reactants["Molar Concentration"] = (reactants["Mole Fraction"]*P_reac)/(Ru*T_reac)

# Fuel to air on mass basis
s = (inv_equiv_ratio*a*O2_kg_kmol + inv_equiv_ratio*a*N2_mole_per_O2*N2_kg_kmol)/fuel_kg_kmol
print(f"{s = }")
# unburned mixture density
rho_reac = (reactants_molecular_weight*P_reac)/(Ru*T_reac)
print(f"{rho_reac = } kg/m^3")
alpha_air = k_avg/(rho_reac*cp_avg) # m^2/s, thermal diffusivity
print(f"{alpha_air = } m^2/s")

# in the reaction zone
rho_avg_reaction = (reactants_molecular_weight*P_reac)/(Ru*T_avg_reaction)
print(f"{rho_avg_reaction = } kg/m^3")
Y_fuel_avg_reaction = 0.5*get_value(reactants, "CxHyOz", "Mass Fraction")
Y_O2_avg_reaction = 0.5*get_value(reactants, "O2", "Mass Fraction")
print(f"{Y_fuel_avg_reaction = }")
print(f"{Y_O2_avg_reaction = }")

K = A * T_avg_reaction**b_coeff * np.exp(-TA/T_avg_reaction)
print(f"{K = }")
fuel_k = -K * rho_avg_reaction**(m+n) \
        * (Y_fuel_avg_reaction/fuel_kg_kmol)**m \
        * (Y_O2_avg_reaction/O2_kg_kmol)**n
print(f"d[CxHyOz]/dt = {fuel_k} kmol/m3-s")
fuel_flux = fuel_k*fuel_kg_kmol
print(f"{fuel_flux = } kg/m3-s")

# Laminar flame speed
sL = np.sqrt(-2*alpha_air*(s+1)*(fuel_flux/rho_reac))
print(f"{sL = } m/s")
# Laminar flame thickness
delta = (2*alpha_air)/sL
print(f"{delta*1000 = } mm")

print(reactants)

print("GLOBAL ONE-STEP APPROXIMATIONS PRODCUTS.")
products = {
    "Gas": ["CO2", "H2O", "O2", "N2"],
    "kmol": [b, c, (inv_equiv_ratio-1)*e, inv_equiv_ratio*f],
    "kg/kmol": [CO2_kg_kmol, H2O_kg_kmol, O2_kg_kmol, N2_kg_kmol],
}
products = pd.DataFrame(products)
products_kmols = products["kmol"].sum()
print(f"Poducts Total kmol: {products_kmols} kmol")
products["Mole Fraction"] = products["kmol"]/products_kmols
products_molecular_weight = (products["Mole Fraction"]*products["kg/kmol"]).sum()
print(f"Poducts Total Molecular Weight: {products_molecular_weight} kg/kmol")
print(f"Poducts Total Mass: {products_kmols*products_molecular_weight} kg")
products["Mass Fraction"] = (products["Mole Fraction"]*products["kg/kmol"])/products_molecular_weight

print("GLOBAL ONE-STEP APPROXIMATIONS PRODCUTS.")
print(products)

