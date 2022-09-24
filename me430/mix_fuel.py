import numpy as np

# MOLE FRACTION of each fuel
x_fuel = np.array([
    0.4,
    0.6,
])

# CxHyOz for each fuel molecule
formula = np.array([
    [7, 16, 0],  # C7H16
    [8, 18, 0],  # C8H18
])
weighted_formula = np.einsum('ij,i->ij',formula,x_fuel)
print(weighted_formula)
print(weighted_formula.sum(0))


