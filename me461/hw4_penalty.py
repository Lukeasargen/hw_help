import numpy as np


a = 1e6 * np.array(
    [[2.5e4, -1, 0],
     [1, 1.5, -0.5],
     [0, -0.5, 0.5]]
)

# Dependent variables
b = np.array(
    [1,
    5,
    2]
)

q = np.linalg.solve(a, b) / 1e-6

print(f"Q0: {q[0]:.3f}")
print(f"Q1: {q[1]:.3f}")
print(f"Q2: {q[2]:.3f}")


