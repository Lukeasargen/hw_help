import numpy as np

# Coefficient Matrix

# ax, ay, bx, by, ff, nn

a = np.array(
    [
        [-1, 1, -1, 1, 0, 0],
        [1, 0, -1, 0, 0, 0],
        [0, 1, 0, -1, 0, 0],
        [0, 0, 1, 0, -1, 0],
        [0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, -1, 0]
    ]
)

# Dependent variables
b = np.array(
    [
        0,
        0,
        1.5*9.81,
        0,
        2.5*9.81,
        (0.2486*-25)/0.33
    ]
)

x = np.linalg.solve(a, b)

print("ax :", x[0])
print("ay :", x[1])
print("bx :", x[2])
print("by :", x[3])
print("ff :", x[4])
print("nn :", x[5])

print("uk :", x[4]/x[5])
