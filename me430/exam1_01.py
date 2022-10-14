
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
    "N2": [0],
    "NO": [90297],
    "NO2": [33098],
    "N2O": [82574],
    "NH": [358790],
}

df = pd.DataFrame(data)
print(df)

# Reaction 1: O + N2O <-> NO + NO
reac = df["NO"] + df["NO"] - df["O"] - df["N2O"]
print(f"R1 Enthalpy of Reaction Forward:\n{reac}")

# Reaction 1: H + N2O <-> NO + NH
reac = df["NO"] + df["NH"] - df["H"] - df["N2O"]
print(f"R2 Enthalpy of Reaction Forward:\n{reac}")

# Reaction 1: O + N2 + M <-> N2O + M
reac = df["N2O"] - df["N2"] - df["O"]
print(f"R3 Enthalpy of Reaction Forward:\n{reac}")
