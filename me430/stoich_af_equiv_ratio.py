
C_kg_kmol = 12
H_kg_kmol = 1
O_kg_kmol = 16
N_kg_kmol = 14

# Fuel
# 1 molar unit
# Cx Hy Oz
x = 2
y = 6
z = 1
equiv_ratio = 0.9

O2_vol = 0.21
N2_vol = 0.79  # includes the other intert gases like Argon

# Calculations below
N2_mole_per_O2 = 3.76  # N2_vol/O2_vol 3.76 moles N2 per 1 mole O2
print(f"{N2_mole_per_O2=}")

N2_kg_kmol = 2*N_kg_kmol
O2_kg_kmol = 2*O_kg_kmol

# CxHyOz + a*(O2+3.76N2) -> b*CO2 + c*H2O + d*N2
a = x + y/4 - z/2
b = x
c = y/2
d = N2_mole_per_O2 * a
print(f"{a=:.4e}")
print(f"{b=:.4e}")
print(f"{c=:.4e}")
print(f"{d=:.4e}")

print("WITHOUT equiv_ratio")
air_kg_kmol = a * (O2_kg_kmol + N2_mole_per_O2 * N2_kg_kmol)
fuel_kg_kmol = x*C_kg_kmol + y*H_kg_kmol + z*O_kg_kmol
AF_stoich = air_kg_kmol/fuel_kg_kmol
print(f"{air_kg_kmol=}")
print(f"{fuel_kg_kmol=}")
print(f"{AF_stoich=}")

print("WITH equiv_ratio")
air_kg_kmol = (1/equiv_ratio) * a * (O2_kg_kmol + N2_mole_per_O2 * N2_kg_kmol)
fuel_kg_kmol = x*C_kg_kmol + y*H_kg_kmol + z*O_kg_kmol
AF_stoich = air_kg_kmol/fuel_kg_kmol
print(f"{air_kg_kmol=}")
print(f"{fuel_kg_kmol=}")
print(f"{AF_stoich=}")


