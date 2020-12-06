import argparse
import numpy as np

"""
Coordinate System

Thin wall pressure vessel.
Tubular shaft.

Positive moment makes the beam smile :)
It is about the Z axis.

            | torque
            v
moment >---------B-----------<
      ( |                   | ) moment
axial   |                   |
force ->|                   | <- axial force
        |                   |
        |                   |
        ---------------------
                          ^
                          | torque
   ^ y+
   |
z+ .--> x+

Location of the element is B is on outer radius.
B is the on the top. Largest y distance.

Take element looking down the y- axis

                | reaction torque
    reaction    v
    moment >---------------------<
          ( |                   | ) moment
    reaction|                   |
    force ->|B                  | <- axial force
            |                   |
            |                   |
            ---------------------
                            ^
                            | torque

.-> x+
|
v z+
"""

parser = argparse.ArgumentParser(description="Shear Tranforms and Mohr's Circle")
parser.add_argument('--ri', type=float, required=True, metavar="Inner Radius")
parser.add_argument('--d', type=float, required=True, metavar="Thickness")
parser.add_argument('--f', type=float, default=0.0, metavar="Axial Force")
parser.add_argument('--t', type=float, default=0.0, metavar="Torque")
parser.add_argument('--m', type=float, default=0.0, metavar="Moment")
parser.add_argument('--p', type=float, default=0.0, metavar="Pressure")
parser.add_argument('--v', action='store_true')
args = parser.parse_args()

inner_radius = args.ri
thickness = args.d

if thickness == 0:
    print("Thickness cannot be 0.")
    exit()

axial_force = args.f
torque = args.t
moment = args.m
pressure = args.p

visualize = args.v

print("Inner Radius = {}".format(inner_radius))
print("Thickness = {}".format(thickness))

print("Axial Force = {}".format(axial_force))
print("Torque = {}".format(torque))
print("Moment = {}".format(moment))
print("Pressure = {}\n".format(pressure))

ratio = inner_radius/thickness
print("Radius/Thickness = {}".format(ratio))
if ratio > 10:
    print("Ratio is >10. Valid thin wall p.v. approximation.")
else:
    print("Ratio is <10. Does not satisfy thin wall p.v. approximations")

# Calculations
outer_radius = inner_radius+thickness
print("\nOuter Radius = {}".format(outer_radius))

area = np.pi*(outer_radius**2 - inner_radius**2)
print("Area = {}".format(area))

J = (np.pi/2) * (outer_radius**4 - inner_radius**4)
print("J = {}".format(J))

I = (np.pi/4) * (outer_radius**4 - inner_radius**4)
print("I = {}\n".format(I))


# Pressure Stresses
hoop_stress = (pressure*inner_radius) / (thickness)
print("Hoop Pressure Stress = {}".format(hoop_stress))
longitudinal_stress = (pressure*inner_radius) / (2*thickness)
print("Longitudinal Pressure Stress = {}".format(longitudinal_stress))

# Axial Stresses
axial_stress = axial_force / area
print("Axial Stress = {}".format(axial_stress))

# Troque stresses
shear = (torque * outer_radius) / J
print("Shear Stress = {}".format(shear))

# Moment / bendig stresses
# negative for a positive moment and positive y in the axis defined above
bending_stress = (-moment*outer_radius) / I
print("Bending Stress = {}".format(bending_stress))


# Sumatation
sigma_axial = longitudinal_stress + axial_stress + bending_stress
sigma_hoop = hoop_stress

print("Total Stresses.")
print("σ axial={:.4f}, σ hoop={:.4f}, τ={:.4f}".format(sigma_axial, sigma_hoop, shear))


if visualize:
    from PIL import Image, ImageDraw, ImageFont

    # Size and position
    size = (1000, 1000)
    scale = 90
    font_size = int(scale / 4.8)
    # print("font_size :", font_size)
    font_path = "fonts/InputSans-Regular.ttf"

    input_center = (0.25*size[0], 0.5*size[1])
    element_center = (0.75*size[0], 0.5*size[1])

    # Colors
    bkg_color = (192, 192, 192)
    axial_color = (222,184,135)
    axial_outline = (0, 0, 0)

    torque_color = (184,222,135)
    torque_outline = (0, 0, 0)

    moment_color = (184,135,222)
    moment_outline = (0, 0, 0)

    vessel_color = (128, 128, 128)
    vessel_outline = (0, 0, 0)
    vessel_width = 4
    text_color = (0,0,0)


    im = Image.new('RGB', size, bkg_color)
    draw = ImageDraw.Draw(im)
    font = ImageFont.truetype(font_path, font_size)

    def square_points(s):
        return np.array([[-s, -s],
                    [ s, -s],
                    [ s,  s],
                    [-s,  s]])

    def arrow_points(s):
        lr = 0.75  # length ratio
        wr = 0.06  # width ratio
        ar = 0.16  # arrow width
        return np.array([[s, 0],  # top
                    [s*lr, s*ar],  # top right out
                    [s*lr, s*wr],  # top right in
                    [-s, s*wr],  # back right in
                    [-s, -s*wr],  # back left in
                    [s*lr, -s*wr],  # top left in
                    [s*lr, -s*ar]])  # top left out

    def rotate_2d(vectors, degrees, center=(0,0)):
        angle = np.radians(degrees)
        c, s = np.cos(angle), np.sin(angle)
        R = np.array(((c, -s), (s, c)))
        out = np.dot(vectors, R) + center
        return out

    def draw_element(center, scale, title="", angle=0, sx=0, sy=0, txy=0):
        # square element
        input_square = rotate_2d(square_points(scale), angle, center)
        draw.polygon(xy=tuple([tuple(x) for x in input_square]), fill=vessel_color, outline=vessel_outline)

        text_move = 1.75

        # shear arrows
        if txy != 0:
            arr = rotate_2d(arrow_points(scale), 90+angle) + rotate_2d([scale*1.25, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            
            arr = rotate_2d(arrow_points(scale), 90+angle) + rotate_2d([scale*1.25, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            arr = rotate_2d(arrow_points(scale), -90+angle) + rotate_2d([-scale*1.25, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            
            arr = rotate_2d(arrow_points(scale), -90+angle) + rotate_2d([-scale*1.25, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            arr = rotate_2d(arrow_points(scale), 0+angle) + rotate_2d([scale*1.25, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            
            arr = rotate_2d(arrow_points(scale), 0+angle) + rotate_2d([scale*1.25, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            arr = rotate_2d(arrow_points(scale), 180+angle) + rotate_2d([-scale*1.25, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)

            arr = rotate_2d(arrow_points(scale), 180+angle) + rotate_2d([-scale*1.25, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)

            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], 225+angle, center)), text="txy={:.2f}".format(txy), fill=text_color, anchor="mm", font=font)

        # normal arrows
        if sx != 0:
            arr = rotate_2d(arrow_points(scale*0.6), 0+angle) + rotate_2d([scale*1.6, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            
            arr = rotate_2d(arrow_points(scale*0.6), 0+angle) + rotate_2d([scale*1.6, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            arr = rotate_2d(arrow_points(scale*0.6), 180+angle) + rotate_2d([-scale*1.6, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            
            arr = rotate_2d(arrow_points(scale*0.6), 180+angle) + rotate_2d([-scale*1.6, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], 180+angle, center)), text="sx={:.2f}".format(sx), fill=text_color, anchor="mm", font=font)

        if sy != 0:
            arr = rotate_2d(arrow_points(scale*0.6), 90+angle) + rotate_2d([scale*1.6, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            
            arr = rotate_2d(arrow_points(scale*0.6), 90+angle) + rotate_2d([scale*1.6, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            arr = rotate_2d(arrow_points(scale*0.6), -90+angle) + rotate_2d([-scale*1.6, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            
            arr = rotate_2d(arrow_points(scale*0.6), -90+angle) + rotate_2d([-scale*1.6, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], -90+angle, center)), text="sy={:.2f}".format(sy), fill=text_color, anchor="mm", font=font)

        # Title
        draw.text(xy=tuple(rotate_2d([-scale*2.25, -scale*2.25], 0, center)), text="{}, θ={:.2f}.".format(title, angle), fill=text_color, anchor="mm", font=font)

    def draw_vessel(center, scale, title="", radius=0, thickness=0, axial_force=0, torque=0, moment=0, pressure=0):
        c = center
        s = scale
        # Draw the cylinder
        r = 2  # roundness
        e = 1 + (r/4)
        draw.line(xy=[c[0]-e*s, c[1]-s, c[0]+e*s, c[1]-s], fill=vessel_color, width=vessel_width)
        draw.line(xy=[c[0]-e*s, c[1]+s, c[0]+e*s, c[1]+s], fill=vessel_color, width=vessel_width)
        draw.arc(xy=[c[0]-r*s, c[1]-s, c[0]-s, c[1]+s], start=90, end=-90, fill=vessel_color, width=vessel_width)
        draw.arc(xy=[c[0]+s, c[1]-s, c[0]+r*s, c[1]+s], start=-90, end=90, fill=vessel_color, width=vessel_width)
        
        # Draw point of interest and cut plane
        draw.text(xy=(c[0], c[1]-s), text="B", fill=text_color, anchor="mm", font=font)
        draw.line(xy=[c[0], c[1]-2.5*s, c[0], c[1]+2.5*s], fill=vessel_color, width=vessel_width)
        draw.text(xy=(c[0], c[1]-2.5*s), text="Cut Plane", fill=text_color, anchor="mm", font=font)

        # Draw the pressure
        if pressure != 0:
            draw.text(xy=c, text="P={:.2f}".format(pressure), fill=text_color, anchor="mm", font=font)

        # Draw the axial_force
        if axial_force != 0:
            arr = rotate_2d(arrow_points(s/3), 180, [c[0]+r*s, c[1]])
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            draw.text(xy=(c[0]+s,c[1]-s/5), text="Axial Force={:.2f}".format(axial_force), fill=text_color, anchor="mm", font=font)
            # Reaction at cut plane
            arr = rotate_2d(arrow_points(s/2), 0, [c[0]-r*s, c[1]])
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=axial_color, outline=axial_outline)
            # draw.text(xy=(c[0]+s,c[1]-s/5), text="Axial Force\nNormal Stress.", fill=text_color, anchor="mm", font=font)

        # Draw the torque
        if torque != 0:
            arr = rotate_2d(arrow_points(s/3), 90, [c[0]+s, c[1]+1.2*s])
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=torque_color, outline=torque_outline)
            draw.text(xy=(c[0]+s,c[1]+1.2*s), text="Torque={:.2f}".format(torque), fill=text_color, anchor="mm", font=font)
            # Reaction at cut plane
            arr = rotate_2d(arrow_points(s/3), 270, [c[0], c[1]-1.2*s])
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=torque_color, outline=torque_outline)
            draw.text(xy=(c[0]+s/10,c[1]-1.5*s), text="Reaction Torque\nShear stress along cut.", fill=text_color, anchor="mm", font=font)

        # Draw the moment
        if moment != 0:
            arr = rotate_2d(arrow_points(s/3), 135, [c[0]+1.9*s, c[1]-0.9*s])
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=moment_color, outline=moment_outline)
            draw.text(xy=(c[0]+1.9*s,c[1]-1*s), text="Moment={:.2f}".format(moment), fill=text_color, anchor="mm", font=font)
            # Reaction at cut plane
            arr = rotate_2d(arrow_points(s/3), 45, [c[0]-s/5, c[1]-0.9*s])
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=moment_color, outline=moment_outline)
            draw.text(xy=(c[0]-2.1*s,c[1]-1.4*s), text="Reaction Moment\nBending Stress.", fill=text_color, anchor="mm", font=font)


        # Title
        draw.text(xy=(c[0] -s*2.25, c[1] -s*2.25), text="{}".format(title), fill=text_color, anchor="mm", font=font)


    # Draw the input pressure vessel
    draw_vessel(center=input_center, scale=scale, title="Input",
                radius=inner_radius,
                thickness=thickness,
                axial_force=axial_force,
                torque=torque,
                moment=moment,
                pressure=pressure)

    # Draw the stress element
    draw_element(center=element_center, scale=0.8*scale, title="Stress Element at B", sx=sigma_axial, sy=sigma_hoop, txy=shear)

    im.show()


