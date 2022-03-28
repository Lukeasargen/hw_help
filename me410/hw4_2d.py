import numpy as np

k = 15.1  # W/m*K, conduction heat transfer coefficient
h = 16  # W/m^2*K, convection heat transfer coefficient
T_b = 95  # C, base temperature
T_inf = 25  # C, surrounding temperature
L = 0.18  # m, length
W = 0.01  # m, width
H = 0.002  # m, height
A_c = W*H  # m^2, cross section area
print(f"{A_c=} m^2")
P = 2*W + 2*H  # m, cross section perimeter
print(f"{P=} m^2")
num_nodes = 2
delta_x = L/num_nodes
print(f"{delta_x=} m")
theta_b = T_b-T_inf
print(f"{theta_b=} C")
# Table 3.4, Case A, tip convection
m_sq = (h*P)/(k*A_c)
m = np.sqrt(m_sq)
temp_ratio = lambda x: ( np.cosh(m*(L-x)) + (h/(m*k))*np.sinh(m*(L-x)) ) / ( np.cosh(m*L) + (h/(m*k))*np.sinh(m*L) )
temp = lambda x: theta_b*temp_ratio(x) + T_inf
q_out = np.sqrt(h*P*k*A_c)*theta_b*( np.sinh(m*L) + (h/(m*k))*np.cosh(m*L) ) / ( np.cosh(m*L) + (h/(m*k))*np.sinh(m*L) )
print(f"{temp_ratio(0.09)=}")
print(f"{temp(0.09)=}")
print(f"{temp_ratio(0.18)=}")
print(f"{temp(0.18)=}")
print(f"{q_out=}")



A = -2 - (h*P*delta_x**2)/(k*A_c)
B = -1 - (h*P*delta_x**2)/(k*A_c*2) - (h*delta_x)/(k)
print(f"{A=}")
print(f"{B=}")

theta_2 = (theta_b)/(A*B-1)
theta_1 = (-B*theta_b)/(A*B-1)
print(f"{theta_1=} C")
print(f"{theta_2=} C")

T_1 = theta_1+T_inf
T_2 = theta_2+T_inf
print(f"{T_1=} C")
print(f"{T_2=} C")

q_out_nodes = h*P*delta_x*(0.5*theta_b+theta_1+0.5*theta_2) + h*A_c*theta_2
print(f"{q_out_nodes=} W")
