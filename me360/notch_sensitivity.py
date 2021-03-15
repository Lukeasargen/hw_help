import numpy as np

# 6-10, 588

use_si = False  # True=metric

r = 1  # notch radius, mm or inches
sut = 105  # MPa or kpsi

# from the book, appendix
# starts at A-15-2 on 1806
kt = 1.68
kts = 1.42

q = 0.85
qs = 0.88

# calculate sqrtr_a, Neuber constant (mm)^1/2 or (inch)^1/2
if use_si:
    if sut > 340 and sut < 1500:
        sqrtr_a = 1.24 - 2.25e-3*sut + 1.6e-6*sut**2 - 4.11e-10*sut**3
        sqrtr_a_s = 0.958 - 1.83e-3*sut + 1.43e-6*sut**2 - 4.11e-10*sut**3
    else:
        print("sut is out the range for these equations")
else:
    if sut > 50 and sut < 220:
        sqrtr_a = 0.246 - 3.08e-3*sut + 1.51e-5*sut**2 - 2.67e-8*sut**3
        sqrtr_a_s = 0.190 - 2.51e-3*sut + 1.35e-5*sut**2 - 2.67e-8*sut**3
    else:
        print("sut is out the range for these equations")   

print("sqrtr_a :", sqrtr_a)
print("sqrtr_a_s :", sqrtr_a_s)

# if q=0, no sensitivity kf=1
# if q=1, full sensitivity, kf=kt
if q == None:
    q = 1 / (1 + (sqrtr_a / np.sqrt(r)))
print("q :", q)

if qs == None:
    qs = 1 / (1 + (sqrtr_a_s / np.sqrt(r)))
print("qs :", qs)

kf = 1 + q*(kt-1)
print("kf :", kf)

kfs = 1 + qs*(kts-1)
print("kfs :", kfs)
