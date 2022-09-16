# Constants
C_kg_kmol = 12
H_kg_kmol = 1
O_kg_kmol = 16
N_kg_kmol = 14

# Fuel
# 1 molar unit
# Cx Hy Oz
x = 7
y = 16
z = 0

# Use to get
#   - mole fractions in reactants
#   - mole fractions in products
#   - stoich air-fuel ratio
equiv_ratio = 0.8  # <1.0

# use to find the equiv_ratio need for each product to have this mole fraction
x_prod = 0.05

# Calculations below
inv_equiv_ratio = 1/equiv_ratio

O2_vol = 0.21
N2_vol = 0.79  # includes the other intert gases like Argon
N2_mole_per_O2 = 3.76  # N2_vol/O2_vol 3.76 moles N2 per 1 mole O2

N2_kg_kmol = 2*N_kg_kmol
O2_kg_kmol = 2*O_kg_kmol

print("CxHyOz + a*(O2+3.76N2) -> b*CO2 + c*H2O + e*O2 + f*N2")
a = x + y/4 - z/2
b = x
c = y/2
e = (inv_equiv_ratio-1)*a
f = N2_mole_per_O2*inv_equiv_ratio*a
print(f"{a=:.4e}")
print(f"{b=:.4e}")
print(f"{c=:.4e}")
print(f"{e=:.4e}")
print(f"{f=:.4e}")

print("Mole Fractions using equiv_ratio")
air_kmol = a*(1+N2_mole_per_O2)
reactant_kmol = 1 + air_kmol
print(f"REACTANTS. {reactant_kmol=}")
print(f"X CxHyOz = {1/reactant_kmol:.4f}")
print(f"X O2+{N2_mole_per_O2}N2 = {air_kmol/reactant_kmol:.4f}")

print("AF with equiv_ratio")
air_kg_kmol = inv_equiv_ratio * a * (O2_kg_kmol + N2_mole_per_O2 * N2_kg_kmol)
fuel_kg_kmol = x*C_kg_kmol + y*H_kg_kmol + z*O_kg_kmol
AF_stoich = air_kg_kmol/fuel_kg_kmol
print(f"{air_kg_kmol=}")
print(f"{fuel_kg_kmol=}")
print(f"{AF_stoich=}")

product_kmol = b+c+e+f
print(f"PRODUCTS. {product_kmol=}")
print(f"X CO2 = {b/product_kmol:.4f}")
print(f"X H2O = {c/product_kmol:.4f}")
print(f"X O2 = {e/product_kmol:.4f}")
print(f"X N2 = {f/product_kmol:.4f}")

print(f"PRODUCTS set to {x_prod=}")
equiv_ratio_O2 = (a*(1-x_prod*(1+N2_mole_per_O2)))/(x_prod*(x+y/2-a)+a)
print(f"X O2=x_prod: equiv_ratio={equiv_ratio_O2}")
