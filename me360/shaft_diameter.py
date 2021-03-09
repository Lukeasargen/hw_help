
import numpy as np

n = 0
d = 50e-3
se = 123.985e6
sut = 700e6
kf = 1
kfs = 1
ma = 800
ta = 0
mm = 0
tm = 600

A = np.sqrt( 4*((kf*ma)**2) + 3*((kfs*ta)**2) )
print("A :", A)

B = np.sqrt( 4*((kf*mm)**2) + 3*((kfs*tm)**2) )
print("B :", B)

if n != 0:
    d = np.cbrt((16*n/np.pi)*((A/se)+(B/sut)))
    print("d :", d)

if d != 0:
    n = (np.pi*(d**3)/16)*(1/((A/se)+(B/sut)))
    print("n :", n, 1/n)
