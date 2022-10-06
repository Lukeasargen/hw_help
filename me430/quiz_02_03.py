
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = {
    "Attibute": ["H"],
    "Temperature": [0],
    "CO": [-110541],
    "OH": [38985],
    "CO2": [-393546,],
    "H": [217977],
    "O2": [0],
    "O": [249197],
    "H2Og": [-241845],
    "H2": [0],
    "NO": [90297],
    "NO2": [33098],
}

df = pd.DataFrame(data)
print(df)

# Reaction 1: H2Og + CO <-> H2 + CO2
reac = df["H2"] + df["CO2"] - df["H2Og"] - df["CO"]
print(f"Reaction 1 Enthalpy:\n{reac}")

# Reaction 2: NO + O <-> NO2
reac = df["NO2"] - df["NO"] - df["O"]
print(f"Reaction 2 Enthalpy:\n{reac}")

# # Reaction 3: CO2 + H <-> CO + OH
reac = df["CO"] + df["OH"] - df["CO2"] - df["H"]
print(f"Reaction 3 Enthalpy:\n{reac}")
