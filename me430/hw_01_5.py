
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = {
    "Gas": ["N2", "CO2", "H2O", "O2", "CO", "NO"],
    "kmol": [10, 1.5, 2, 0.5, 0.1, 0.005],
    "kg/kmol": [14, 44, 18, 32, 28, 30],
}

df = pd.DataFrame(data)

total_kmols = df["kmol"].sum()
print(f"Total kmol: {total_kmols}")

df["Mole Fraction"] = df["kmol"]/total_kmols
total_molecular_weight = (df["Mole Fraction"]*df["kg/kmol"]).sum()
print(f"Total Molecular Weight: {total_molecular_weight}")

print(f"Total Mass: {total_kmols*total_molecular_weight}")

df["Mass Fraction"] = (df["Mole Fraction"]*df["kg/kmol"])/total_molecular_weight
print(df)

