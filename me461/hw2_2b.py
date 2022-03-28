import matplotlib.pyplot as plt
import numpy as np

num_points = 100
L = 1
E = 1
I = 1
P = 1

xdata = np.linspace(0,L, num_points)

exact = ((P*L)/(6*E*I))*(3-(xdata/L))*xdata**3
ritz = (P/(E*I))*(((L*xdata**2)/2)-((xdata**3)/6))

fig = plt.figure(figsize=(8,4))
ax = fig.subplots(nrows=1, ncols=1)
ax.plot(xdata, exact, label="Exact")
ax.plot(xdata, ritz, label="Ritz Method")
ax.legend()
plt.show()