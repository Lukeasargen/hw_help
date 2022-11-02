# water-gas shift
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Constants
C_kg_kmol = 12
H_kg_kmol = 1
O_kg_kmol = 16
N_kg_kmol = 14
Ru = 8315  # J/Kmol*K, reactants

# Inputs
# Fuel, 1 molar unit, CxHyOz
x = 3; y = 8; z = 0
equiv_ratio = 1.2  # <1.0
P = 4*101325  # atm -> Pa, reactants
T_reac = 1000  # K, reactants
# water-gas shift reactions
T_wgsr = 3000  # K
Kp = 0.138  # parital-pressure equilibrium constant at the given temperature

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
CO_kg_kmol = C_kg_kmol + O_kg_kmol
H2O_kg_kmol = 2*H_kg_kmol + O_kg_kmol
H2_kg_kmol = 2*H_kg_kmol

# Global one-step mechanism, assume no O2 in products
print("CxHyOz + a*(O2+3.76N2) -> b*CO2 + c*CO + d*H2O + e*H2 + f*N2")
print(f"Water-Gas Shift Reaction: Kp(T={T_wgsr:.0f}K) = {Kp}")
# 1) Find a from LHS
a = inv_equiv_ratio*(x + y/4 - z/2)
print(f"{a=}")

# 2) Find Kp for T, this is given

# 3) Solve for b, 0<b<x, Eq 2.74 from Turns
part1 = (((2*a+z)*(Kp-1))+x+(y/2)) / (2*(Kp-1))
part2 = ((((2*a+z)*(Kp-1))+x+(y/2))**2 - (4*Kp*(Kp-1)*(x*z + 2*a*x - x**2)))**0.5 / (2*(Kp-1))
b = part1-part2  # Should be a +-, but the negative root is correct
print(f"{b=}")

# 4) Solve for c, d, e, f
c = x - b
d = 2*a + z - b - x
e = y/2 - d
f = 3.76*a
print(f"{c=}")
print(f"{d=}")
print(f"{e=}")
print(f"{f=}")

reactants = {
    "Gas": ["CxHyOz", "O2", "N2"],
    "kmol": [1, a, a*N2_mole_per_O2],
    "kg/kmol": [fuel_kg_kmol, O2_kg_kmol, N2_kg_kmol],
}
reactants = pd.DataFrame(reactants)

products = {
    "Gas": ["CO2", "CO", "H2O", "H2", "N2"],
    "kmol": [b, c, d, e, f],
    "kg/kmol": [CO2_kg_kmol, CO_kg_kmol, H2O_kg_kmol, H2_kg_kmol, N2_kg_kmol],
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
reactants["Molar Concentration"] = (reactants["Mole Fraction"]*P)/(Ru*T_reac)
print(reactants)

print("GLOBAL FUEL-RICH ONE-STEP APPROXIMATIONS PRODCUTS.")
print(products)
