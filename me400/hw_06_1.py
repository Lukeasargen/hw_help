# gas turbine, 2 comp, 2 turb, reheat, intercooling, regen
# air
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

pressure_ratio = 4  # P2/P1 = P4/P3
q56 = q78 = 300*1000  # J/kg heat add in combustion
regen_add_temp = 20  # C, after second comp

P1 = 100  # kPa
T1 = 290.15  # K

T2 = T1*(pressure_ratio)**((gamma-1)/gamma)
P2 = pressure_ratio*P1
print(f"{T2=} K")
print(f"{P2=} kPa")

T3 = T1
P3 = P2
print(f"{T3=} K")
print(f"{P3=} kPa")

T4 = T2
P4 = pressure_ratio*P3
print(f"{T4=} K")
print(f"{P4=} kPa")

T5 = T4 + regen_add_temp
P5 = P4
print(f"{T5=} K")
print(f"{P5=} kPa")

T6 = q56/cp + T5
P6 = P4
print(f"{T6=} K")
print(f"{P6=} kPa")

T7 = T5
P7 = P6*(T7/T6)**(gamma/(gamma-1))
print(f"{T7=} K")
print(f"{P7=} kPa")

T8 = T6
P8 = P7
print(f"{T8=} K")
print(f"{P8=} kPa")

T9 = T5
P9 = P1
print(f"{T9=} K")
print(f"{P9=} kPa")

eta_th = ((T6-T7)+(T8-T9)-(T2-T1)-(T4-T3))/((T6-T5)+(T8-T7))
print(f"{eta_th=}")
