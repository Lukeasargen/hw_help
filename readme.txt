hookes_law.py
    --v Possion Ratio.
    --e Elastic Modulus.
    --sx Stress X Axis.
    --sy Stress Y Axis.
    --sz Stress Z Axis.
    Example: python hookes_law.py --v 0.32 --200e9 --sx 32000 --sy -10000 --sz -15000

pressure_vessel.py
    --ri Inner Radius.
    --d Wall Thickness.
    --f Axial Force. Positive Compression
    --t Torque. Use --v to check direction.
    --m Moment. Use --v to check direction.
    --p Pressure.
    --v Generate Visual.
    Example: python pressure_vessel.py --ri 1.25 --d 0.02 --f 100 --t -50 --m 35 --p 10 --v

stress_element_mohrs.py
    --sx Stress X Axis.
    --sy Stress Y Axis.
    --txy Shear Stress.
    --deg Angle of cut plane positive in ccw direction from X axis in DEGREES. This output is only printed and is not shown in the visualization.
    --v Generate Visual.
    Example: python stress_element_mohrs.py --sx 125 --sy -13 --txy 44 --deg 15 --v

