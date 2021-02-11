from math import sqrt, pi, degrees, acos

# cos phi = r / z
# xl = omgea L
# xc = 1 / omega C
# omega = 2 pi hx = 2 pi / T

resistor = 2.1
inductor = 2e-3
capacitor = 2.5e-3
hz = 0
omega = 1000
# hz = 1 / T
xc = 0
xl = 0

if hz != 0:
    xl = 2 * pi * hz * inductor
    if capacitor != 0: # divide by zero
        xc = 1 / (2 * pi * hz * capacitor)

if omega != 0:
    xl = omega * inductor
    if capacitor != 0:  # divide by zero
        xc = 1 / (omega * capacitor)

impedance = sqrt((resistor * resistor) + ((xl-xc)*(xl-xc)))
if impedance == 0:
    print("Impedance is zero????")
    exit()

pf = resistor / impedance
angle = degrees(acos(pf))
print("Impedance :", impedance)
print("Power factor :", pf)
print("Phase Angle Degrees :", angle)
