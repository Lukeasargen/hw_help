import numpy as np

L =  1e-3 * np.array([800, 600, 400])  # m, length
A =  1e-6 * np.array([2400, 1200, 600])  # m^2, area
E =  1e9 * np.array([83, 70, 200])  # Pa, elastic modulus
alpha = 1e-6 * np.array([18.9, 23, 11.7])  # 1/C
delta_temp = 80
P1 = 60e3  # N
P2 = 75e3  # N

print(f"{L=}")
print(f"{A=}")
print(f"{E=}")
print(f"{alpha=}")
print(f"{P1=}")
print(f"{P2=}")

K = E*A/L
print(f"{K*1e-6=}")

theta = E*A*alpha*delta_temp
print(f"{theta=}")

k1, k2, k3 = K
K = np.array(
    [[k1, -k1, 0, 0],
     [-k1, k1+k2, -k2, 0],
     [0, -k2, k2+k3, -k3],
     [0, 0, -k3, k3]]
)
print(f"{K*1e-6=}")

t1, t2, t3 = theta

F = np.array(
    [-t1,
    t1-t2-P1,
    t2-t3-P2,
    t3]
)
print(f"{F*1e-6=}")

q2, q3 = np.linalg.solve(K[1:3, 1:3], F[1:3])
q = np.array([0, q2, q3, 0])
print(f"{q*1e-3=}")


for i in range(3):
    stress = E[i] * ( (-q[i]+q[i+1])/L[i] - (alpha[i]*delta_temp) )
    print(f"Element {i+1}: {stress/1e6} MPa")
    

print(f"R1 = {np.dot(K[0], q)-F[0]}")
print(f"R4 = {np.dot(K[3], q)-F[3]}")


