import numpy as np
import matplotlib.pyplot as plt


x_min = 0
x_max = 1e6

x = np.linspace(x_min, x_max, int(x_max/10))

y1 = 20*np.log10(1/np.sqrt(1+(x/100)**2))
y2 = 20*np.log10(1/np.sqrt(1+(x/1000)**2))
y3 = 20*np.log10(1/np.sqrt(1+(x/10000)**2))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot(x, y1, label="100 Hz")
ax.plot(x, y2, label="1000 Hz")
ax.plot(x, y3, label="10000 Hz")

ax.set_xscale('log')
ax.set_ylim(-30, 1)
ax.set_xlim(x_min, x_max)

ax.set_ylabel("Magnitude (dB)")
ax.set_xlabel("Frequency (Hz)")

ax.grid()

ax.legend()

plt.show()