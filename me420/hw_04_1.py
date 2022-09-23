import matplotlib.pyplot as plt
import numpy as np


mach_number = np.linspace(1, 10, 100)
print(mach_number)
mach_angle = np.arcsin(1/mach_number)
print(mach_angle)

fig, ax = plt.subplots(figsize=(4, 3))
ax.set_title("Mach Angle as a funciton of Mach Number",fontsize=10)
ax.set_xlabel("Mach Number")
ax.set_ylabel("Mach Angle [rad]")
ax.set_ylim(bottom=0)
ax.plot(mach_number, mach_angle)
ax.grid()
fig.tight_layout()
plt.show()
