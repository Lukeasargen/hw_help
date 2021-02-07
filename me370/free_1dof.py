import numpy as np
import matplotlib.pyplot as plt

# System
m = 2  # k
k = 8  # N/m, stiffness
c = 0  # Ns/m, damping
print("System: m={} kg. k={} N/m. c={} Ns/m.".format(m, k, c))

x0 = 3e-3 # m, inital displacement
v0 = 8e-3 # m/s, intial velocity
print("Initial conditions: x0={} m. v0={} m/s..".format(x0, v0))

omega_n = np.sqrt(k/m)  # natural frequency
print("Natural Frequency = {:.4f} rad/s.".format(omega_n))
amplitude = np.sqrt(x0**2 + (v0/omega_n)**2)
print("Amplitude = {:.4f} m.".format(amplitude))
phase = np.arctan2(x0, v0/omega_n)
print("Phase = {:.4f} rad.".format(phase))

# Solution
x_t = lambda t: amplitude*np.sin(omega_n*t + phase)

sim_length = 8 # seconds
steps = np.linspace(0, sim_length, sim_length*100)

plt.plot(steps, x_t(steps))

plt.grid()
# plt.grid(b=True, which='major', color='#666666', linestyle='-')

plt.xlabel('Seconds')
plt.ylabel('Displacement (m)')

plt.show()


