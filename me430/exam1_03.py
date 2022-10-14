import numpy as np

# Constants
C_kg_kmol = 12
H_kg_kmol = 1
O_kg_kmol = 16
N_kg_kmol = 14

# Fuel, 1 molar unit, CxHyOz
x = 3
y = 6
z = 0

equiv_ratio = 0.9  # <1.0

# Calculations below
inv_equiv_ratio = 1/equiv_ratio
O2_vol = 0.21
N2_vol = 0.79  # includes the other intert gases like Argon
N2_mole_per_O2 = 3.76  # N2_vol/O2_vol 3.76 moles N2 per 1 mole O2
N2_kg_kmol = 2*N_kg_kmol
O2_kg_kmol = 2*O_kg_kmol

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

print("Mole Fractions using equiv_ratio")
O2_kmol = inv_equiv_ratio*a*1
N2_kmol = inv_equiv_ratio*a*N2_mole_per_O2
fuel_kmol = 1
print(f"{O2_kmol=}")
print(f"{N2_kmol=}")
print(f"{fuel_kmol=}")
reactant_kmol = fuel_kmol + O2_kmol + N2_kmol
print(f"REACTANTS. {reactant_kmol=}")
print(f"X CxHyOz = {fuel_kmol/reactant_kmol:.4f}")
print(f"X O2 = {O2_kmol/reactant_kmol:.4f}")
print(f"X N2 = {N2_kmol/reactant_kmol:.4f}")

product_kmol = b+c+e+f
print(f"PRODUCTS. {product_kmol=}")
print(f"X CO2 = {b/product_kmol}")
print(f"X H2O = {c/product_kmol}")
print(f"X O2 = {e/product_kmol}")
print(f"X N2 = {f/product_kmol}")

# Find the molar concentration using ideal gas
P = 5*101325  # 5 atm -> Pa
Ru = 8315  # J/Kmol*K
T = 1200  # K
fuel_c = (fuel_kmol*P)/(reactant_kmol*Ru*T)
O2_c = (O2_kmol*P)/(reactant_kmol*Ru*T)
N2_c = (N2_kmol*P)/(reactant_kmol*Ru*T)
print("MOLAR CONCENTRATION REACTANTS.")
print(f"[CxHyOz] = {fuel_c = } kmol/m3")
print(f"[O2] = {O2_c = } kmol/m3")
print(f"[N2] = {N2_c = } kmol/m3")

# Find the intitial-time-rate-of-change of fuel molar concentration
print("initial time-rate-of-change of fuel molar concentration")
# Values are given, they also come from Table 5.1
A = 2.36e9  # (kmol/m3)^-0.75 / s 
m = -0.1
n = 1.85
TA = 15098  # K

# negative because it's depleted
fuel_k = -A * np.exp(-TA/T) * fuel_c**m * O2_c**n
print(f"d[CxHyOz]/dt = {fuel_k = } kmol/m3-s")


