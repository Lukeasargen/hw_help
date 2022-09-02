
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = {
    "Gas": ["Nitrogen", "Carbon Dioxide", "Helium", "Hydrogen", "Methane"],
    "Heat Capacity Ratio": [1.40, 1.28, 1.66, 1.41, 1.32],  # https://www.engineeringtoolbox.com/specific-heat-ratio-d_608.html
    "Gas Constant [J/kg*K]": [296.8, 188.92, 2077.1, 4124.2, 518.28],  # https://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html
    "c_p(T=300K) [J/kg*K]": [1039, 846, 5192.6, 14307, 2253.7],
    "c_v(T=300K) [J/kg*K]": [743, 657, 3115.6, 10183, 1735.4],
}

df = pd.DataFrame(data)

temperature_c = np.linspace(-10, 40, 6)
temperature_k = temperature_c + 273.15
print(temperature_c)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title("Speed of Sound vs. Temperature")
ax.set_xlabel("Temperature [C]")
ax.set_ylabel("Speed of Sound [m/s]")

for idx, row in df.iterrows():
    # Ideal gas speed of sound
    speed = np.sqrt(temperature_k*row["Heat Capacity Ratio"]*row["Gas Constant [J/kg*K]"])
    # print(speed)
    ax.plot(temperature_c, speed,
        label=f"{row['Gas']} ($\gamma$={row['Heat Capacity Ratio']}, R={row['Gas Constant [J/kg*K]']})")

ax.grid()
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
fig.tight_layout()
plt.show()
