# gas turbine, 2 comp, 2 turb, reheat, intercooling, regen
# air
gamma = 1.4
R = 287.0  # J/Kg*K
cp = 1005  # J/Kg*K

pressure_ratio = 3  # P2/P1 = P4/P3 = P6/P7 = P8/P9
T1 = T3 = 300  # K, compressor enterance
T6 = T8 = 1200  # K, turbine enterance

T2 = T4 = 410.1  # A17
print(f"{T2=} K")
print(f"{T4=} K")

T7 = T9 = 911.95  # A17
print(f"{T7=} K")
print(f"{T9=} K")

eta_th = ((T6-T7)+(T8-T9)-(T2-T1)-(T4-T3))/((T6-T4)+(T8-T7))
print(f"2a) {eta_th=}")
eps = (T2-T1)/(T6-T7)
print(f"2a) {eps=}")

eps = 0.75
T5 = T4 + eps*(T7-T4)
print(f"{T5=} K")

eta_th = ((T6-T7)+(T8-T9)-(T2-T1)-(T4-T3))/((T6-T5)+(T8-T7))
print(f"2b) {eta_th=}")

