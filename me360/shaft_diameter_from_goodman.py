
import numpy as np

n = 1.5  # saftey factor
d = 0  # diameter

sy = 82e3  # yielding
se = 13.09e3  # fatigue
sut = 68e3  # ultimate
kf = 1.7  # normal stress-concentration factor
kfs = 1.5  # shear stress-concentration factor
ma = 3651  # alternating moment
ta = 0  # alternating torsion
mm = 0  # mean moment
tm = 3240  # mean torsion

A = np.sqrt( 4*((kf*ma)**2) + 3*((kfs*ta)**2) )
print("A :", A)

B = np.sqrt( 4*((kf*mm)**2) + 3*((kfs*tm)**2) )
print("B :", B)

if n != 0:
    d = np.cbrt((16*n/np.pi)*((A/se)+(B/sut)))
    print("d :", d)
elif d != 0:
    n = (np.pi*(d**3)/16)*(1/((A/se)+(B/sut)))
    print("from class slide 10 n :", n)
    if n < 1:
        print("This part fails. Use Basquin's equation (S=aN^b) to determine part life.")

print("Additional stress calculations and checks.")
vm_a = np.sqrt( (((32*kf*ma)/(np.pi*d**3))**2) + 3*(((16*kfs*ta)/(np.pi*d**3))**2) )
print("vm_a :", vm_a)

vm_m = np.sqrt( (((32*kf*mm)/(np.pi*d**3))**2) + 3*(((16*kfs*tm)/(np.pi*d**3))**2) )
print("vm_m :", vm_m)

nn = (vm_a/se) + (vm_m/sut)
print("goodman 1/n :", nn)
print("goodman n :", 1/nn)

vm_max = np.sqrt((((32*kf*(mm+ma))/(np.pi*(d**3)))**2) + 3*(((16*kfs*(tm+ta))/(np.pi*(d**3)))**2))
print("vm_max :", vm_max)

ny = sy/vm_max
print("ny=sy/vm_max :", ny)
