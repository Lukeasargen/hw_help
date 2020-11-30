Requirements: pip install numpy pillow


Notes:
    - Compression = Positive Stress
    - Use the same units for all inputs


hookes_law.py
    Prints out the strain values for general hooke's law. 
    --v Possion Ratio.
    --e Elastic Modulus.
    --sx Stress X Axis.
    --sy Stress Y Axis.
    --sz Stress Z Axis.
    Example: python hookes_law.py --v 0.32 --e 200e9 --sx 32000 --sy -10000 --sz -15000


pressure_vessel.py
    Calculates the stress for an element on the surface of a pressure vessel. It also can calculate the stress for axial force, moment, and torque applied to a tubular shaft. For solid shafts you can use 0 for inner radius and set the thickness to the radius. Use --v to see the direction conventions.
    --ri Inner Radius.
    --d Wall Thickness.
    --f Axial Force. Positive Compression
    --t Torque. Use --v to check direction.
    --m Moment. Use --v to check direction.
    --p Pressure.
    --v Generate Visual.
    Example: python pressure_vessel.py --ri 1.25 --d 0.02 --f 100 --t -50 --m 35 --p 10 --v


stress_element_mohrs.py
    Calculates the stress tranformation equations for a given stress element. Use the --deg argument to find the stress at a cut plane given the angle of the plane. 
    --sx Stress X Axis.
    --sy Stress Y Axis.
    --txy Shear Stress.
    --deg Angle of cut plane positive in ccw direction from X axis in DEGREES. This output is only printed and is not shown in the visualization.
    --v Generate Visual.
    Example: python stress_element_mohrs.py --sx 125 --sy -13 --txy 44 --deg 15 --v

