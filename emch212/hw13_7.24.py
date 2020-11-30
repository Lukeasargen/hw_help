import numpy as np

# Coefficient Matrix

# P, Nc, Ne, Nf, Dx, Dy

a = np.array(
    [
        [-2, 0, 0, 0, 1, 0],
        [0, 2, 2, 0, 0, 1],
        [-50/12, -96/12, 104/12, 0, 1, 89/12],
        [0, 0, 0, 0, -1, 0],
        [0, 0, 0, 2, 0, 1],
        [0, 0, 0, 0, -28/12, 110/12]


    ]
)

# Dependent variables
b = np.array(
    [
        (3300/32.2)*(-88/18),
        3300,
        0,
        (4300/32.2)*(-88/18),
        4300,
        0
    ]
)

x = np.linalg.solve(a, b)

print("P :", x[0])
print("Nc :", x[1])
print("Ne :", x[2])
print("Nf :", x[3])
print("Dx :", x[4])
print("Dy :", x[5])

