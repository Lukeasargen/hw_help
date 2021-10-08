
import math


measure_hz = 3000
noise_hz = 100

# need a high-pass filter

# first order filter
# attenuate by more than 80% at 100Hz, <0.2
# attenuate by less than 10% at 3000Hz, >0.9

def high_pass_gain(input_hz, cutoff_hz):
    return 1/ math.sqrt(1+(cutoff_hz/input_hz)**2)

cutoff_hz = 1000
measure_gain = high_pass_gain(measure_hz, cutoff_hz)
noise_gain = high_pass_gain(noise_hz, cutoff_hz)

print(f"Gain at {measure_hz=}: {measure_gain>0.9} {measure_gain=}")
print(f"Gain at {noise_hz=}: {noise_gain<0.2} {noise_gain=}")

capictor = 1e-6  # farads
# cutoff_hz = 1 / 2*math.pi*R*C
# r = 1 / 2*pi*C*cutoff_hz

resistor = 1 / (2*math.pi*capictor*cutoff_hz)
print(f"{resistor=}")

phase = 180/math.pi *math.atan(cutoff_hz/measure_hz)
print(f"{phase=}")



