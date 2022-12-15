import cantera as ct

gas = ct.Solution('gri30.yaml')
air = "O2:0.21,N2:0.79"
phi = gas.equivalence_ratio()
print(f"phi = {phi:1.3f}")

Z = gas.mixture_fraction(fuel="H2:1", oxidizer=air, element="H")
print(f"Z(mixture fraction based on H) = {Z:1.3f}")

print(f"mass fraction of H2 = {gas['H2'].Y[0]:1.3f}")
