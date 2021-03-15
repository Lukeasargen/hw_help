import numpy as np


in_units = input("What units? (si or us) : ")

use_si = True if in_units=="si" else False
if use_si:
    print("Using si units (MPa, m, etc.)")
    sut = float(input("Sut in MPa : "))
    se_p = 0.5*sut if sut < 1400 else 700
    print("se'= {} MPa".format(se_p))
else:
    print("Using us units (kpsi, in, etc.)")
    sut = float(input("Sut in kpsi : "))
    se_p = 0.5*sut if sut < 200 else 100
    print("se'= {} kpsi".format(se_p))

# curve fit parameters for surface factor
table_6_2 = [
    [1.21, 1.38, -0.067],
    [2.0, 3.04, -0.217],
    [11.0, 38.6, -0.65],
    [12.7, 54.9, -0.758]
]

finish = int(input("Finish (Options: 1=ground, 2=machined or cold-drawn, 3=hot-rolled, 4=aas-forged) : "))
if finish in range(1,5):
    a = table_6_2[finish-1][1 if use_si else 0]
    b = table_6_2[finish-1][2]
    print("a=", a)
    print("b=", b)
else:
    print("Invalid. Enter number between 1-4.")
    exit()

ka = a*sut**b
print("Surface factor ka=", ka)

print("Size factor. If nonrotating, use table 6-3 for de and find kb with eq 6-19.")
if use_si:
    in_d = float(input("Diameter in mm (enter 0 if unknown) : "))
    if in_d < 7.62:
        kb = 1
    elif in_d < 51:
        kb = (in_d/7.62)**-0.107
    elif in_d < 254:
        kb = 1.51*in_d**-0.157
    else:
        print("Diameter is larger than 254 mm. Cannot use the curve fit equations.")
        exit()
else:
    in_d = float(input("Diameter in inches (enter 0 if unknown) : "))
    if in_d < 0.3:
        kb = 1
    elif in_d < 2:
        kb = (in_d/0.3)**-0.107
    elif in_d < 10:
        kb = 0.91*in_d**-0.157
    else:
        print("Diameter is larger than 10 in. Cannot use the curve fit equations.")
        exit()

print("Size factor kb=", kb)

table_kc = [
    0.59,
    0.85, 
    1.0
]
loading = int(input("Loading (Options: 1=torsion, 2=axial, 3=bending) : "))
kc = table_kc[loading-1]

print("Load factor kc=", kc)

if use_si:
    temp_in = float(input("Temperature in celsius (enter 0 for room temp) : "))
    if temp_in == 0:
        kd = 1
    else:
        kd = 0.99 + 5.9e-4*temp_in - 2.1e-6*temp_in**2
else:
    temp_in = float(input("Temperature in fahrenheit (enter 0 for room temp) : "))
    if temp_in == 0:
        kd = 1
    else:
        kd = 0.99 + 3.5e-4*temp_in - 6.3e-7*temp_in**2

print("Temperature factor kd=", kd)

table_6_4 = {
    "50": 1.0,
    "90": 0.897,
    "95": 0.868,
    "99": 0.814,
    "99.9": 0.753,
    "99.99": 0.702
}
reliability = input("Reliability % (Options: 50, 90, 95, 99, 99.9, 99.99) : ")
ke = table_6_4[reliability]
print("Reliability factor ke=", ke)

se = se_p*ka*kb*kc*kd*ke
if use_si:
    print("Se= {} MPa".format(se))
else:
    print("Se= {} kpsi".format(se))
