import numpy as np

n = 1
sy = 57e3
T = 2819
d = 1
w = 1/4
h = 1/4

# shearing
ls = (2*np.sqrt(3)*T*n)/(w*sy*d)
print("shear length :", ls)

lc = (4*T*n)/(h*sy*d)
print("crushing length :", lc)

