
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = {
    "Attibute": ["H", "G", "G"],
    "Temperature": [0, 1800, 2000],
    "CO": [-110541, -269164, -285948],
    "OH": [38985, 11958, 9194],
    "CO2": [-393546, -396425, -396410],
    "H": [217977, 118830, 106869],
    "O2": [0, 0, 0],
    "O": [249197, 135049, 121709],
    "H2O": [-241845, -147216, -135643]
}

df = pd.DataFrame(data)
print(df)

# Reaction 1: CO + OH <-> CO2 + H
reac = df["CO2"] + df["H"] - df["CO"] - df["OH"]
print(f"Reaction 1 Enthalpy and Gibbs:\n{reac}")
kp = np.exp( -reac[1:]/(8.314*df["Temperature"][1:]) )
print(f"Reaction 1 partial-pressure equilibrium constant:\n{kp}")

# Reaction 2: CO + O2 <-> CO2 + O
reac = df["CO2"] + df["O"] - df["CO"] - df["O2"]
print(f"Reaction 2 Enthalpy and Gibbs:\n{reac}")
kp = np.exp( -reac[1:]/(8.314*df["Temperature"][1:]) )
print(f"Reaction 2 partial-pressure equilibrium constant:\n{kp}")

# Reaction 3: O + H2O <-> OH + OH
reac = df["OH"] + df["OH"] - df["H2O"] - df["O"]
print(f"Reaction 3 Enthalpy and Gibbs:\n{reac}")
kp = np.exp( -reac[1:]/(8.314*df["Temperature"][1:]) )
print(f"Reaction 3 partial-pressure equilibrium constant:\n{kp}")

# Reaction 4: H + O2 <-> OH + O
reac = df["O"] + df["OH"] - df["H"] - df["O2"]
print(f"Reaction 4 Enthalpy and Gibbs:\n{reac}")
kp = np.exp( -reac[1:]/(8.314*df["Temperature"][1:]) )
print(f"Reaction 4 partial-pressure equilibrium constant:\n{kp}")

