import numpy as np

# Coefficient Matrix
# ax, ay, fs
a = np.array(
    [
        [1, 0, 1],
        [0, 1, 0],
        [0, 0, 3],
    ]
)

# Dependent variables
r = 3
l = 4
p = 1000
g = 9.81
fx = p*g*(r/2)*r*l
fy = p*g*r*r*l
w = p*g*(1/4)*np.pi*(r**2)*l
print("fx :", fx)
print("fy :", fy)
print("w :", w)
b = np.array(
    [
        fx,
        w - fy,
        fx*(2*r)/3 + fy*(r/2) - w*(4*r)/(3*np.pi)
    ]
)

x = np.linalg.solve(a, b)
print("ax :", x[0])
print("ay :", x[1])
print("fs :", x[2])
