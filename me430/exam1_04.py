import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Constants
C_kg_kmol = 12
H_kg_kmol = 1
O_kg_kmol = 16
N_kg_kmol = 14

# Fuel, 1 molar unit, CxHyOz
x = 3
y = 8
z = 0

equiv_ratio = 0.95  # <1.0

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
e = (inv_equiv_ratio-1)*a
f = N2_mole_per_O2*inv_equiv_ratio*a
print(f"{a=:.4e}")
print(f"{inv_equiv_ratio*a=:.4e}")
print(f"{b=:.4e}")
print(f"{c=:.4e}")
print(f"{e=:.4e}")
print(f"{f=:.4e}")


reactants = {
    "Gas": ["CxHyOz", "O2", "N2"],
    "kmol": [1, inv_equiv_ratio*a*1, inv_equiv_ratio*a*N2_mole_per_O2],
    "kg/kmol": [fuel_kg_kmol, O2_kg_kmol, N2_kg_kmol],
}
reactants = pd.DataFrame(reactants)

products = {
    "Gas": ["CO2", "H2O", "O2", "N2"],
    "kmol": [b, c, e, f],
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
    print(df)
