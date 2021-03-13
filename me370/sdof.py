import numpy as np

# Equivalent System
m = 1
k = 4
c = 1
print("System: m={}, c={}, k={}".format(m, c, k))

# Initial Conditions
x0 = -0.1
v0 = 0.5
print("Intial: x0={}, v0={}".format(x0, v0))

f0 = 5.0  # forced vibration amplitude
y0 = 0.0  # base excitation amplitude
w = 1.2  # driving frequency

# Calculate system parameters
wn = np.sqrt(k/m)
print("Natural Frequency :", wn)
cc = 2*np.sqrt(k*m)
print("Critical Damping :", cc)
zeta = c/cc
print("Damping ratio :", zeta)
if zeta < 1:
    wd = wn*np.sqrt(1-zeta**2)
    print("Dampened Natural Frequency :", wd)

r = w/wn  # frequency ratio
print("Frequency ratio :", r)

# Generate the response equations
if f0 != 0:  # Forced vibration

    print("EOM: mx'' + cx' + kx = amp * sin(w*t)")
    print("Forcing: f(t) = {} * sin({} * t)".format(f0, w))

    if zeta == 0:
        print("FIXED-BASE EXTERNAL HARMONIC EXCITATION: Undamped")
        if w != wn:
            print("Not at resonsance. w != wn.")
            a2 = f0/(k*(1-r**2))
            print("x(t) = {:.06f} * cos({:.06f} * t) + {:.06f} * sin({:.06f} * t)".format(x0, (v0-a2*w)/wn, wn, a2, w))
        else:
            print("At resonace. w = wn.")
            b2 = (f0*wn)/(2*k)
            print("x(t) = {:.06f} * cos({:.06f} * t) + {:.06f} * t * sin({:.06f} * t)".format(x0, (v0+b2)/wn, wn, b2, w))

    elif zeta < 1:
        print("FIXED-BASE EXTERNAL HARMONIC EXCITATION: Underdamped")
    
    elif zeta == 1:
        print("FIXED-BASE EXTERNAL HARMONIC EXCITATION: Critically Damped")
    
    elif zeta > 1:
        print("FIXED-BASE EXTERNAL HARMONIC EXCITATION: Overdamped")



elif y0 != 0:  # base excitation
    pass
else: # free vibration

    print("EOM: mx'' + cx' + kx = 0")

    if zeta == 0:
        print("FREE VIBRATION: Undamped")
        A = np.sqrt((x0)**2 + (v0/wn)**2)
        phi = np.arctan2(x0*wn, v0)
        print("x(t) = {:.06f} * sin({:.06f} * t + {:.06f})".format(A, wn, phi))

    elif zeta < 1:
        print("FREE VIBRATION: Underdamped")
        A = np.sqrt((x0)**2 + ((v0+zeta*wn*x0)/wd)**2)
        phi = np.arctan2(x0*wd, v0+zeta*wn*x0)
        exp = zeta*wn
        print("x(t) = {:.06f} * e^(-{:.06f} * t) * sin({:.06f}*t + {:.06f})".format(A, exp, wd, phi))

    elif zeta == 1:
        print("FREE VIBRATION: Critically Damped")
        print("x(t) = {:.06f} * e^({:.06f} * t) + {:.06f} * t * e^({:.06f} * t)".format(x0, -zeta*wn, v0+x0*wn, -wn))

    elif zeta > 1:
        print("FREE VIBRATION: Overdamped")
        s1 = (-zeta + np.sqrt(zeta**2 - 1))*wn
        s2 = (-zeta - np.sqrt(zeta**2 - 1))*wn
        a1 = (x0*s2 - v0)/(s2-s1)
        a2 = (v0 - x0*s1)/(s2-s1)
        print("x(t) = {:.06f} * e^({:.06f} * t) + {:.06f} * e^({:.06f} * t)".format(a1, s1, a2, s2))
