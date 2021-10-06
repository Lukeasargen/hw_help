import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Shear Tranforms and Mohr's Circle")
parser.add_argument('--sx', type=float, default=0.0, help="Stress in X axis")
parser.add_argument('--sy', type=float, default=0.0, help="Stress in Y axis")
parser.add_argument('--txy', type=float, default=0.0, help="Shear Stress.")
parser.add_argument('--deg', type=float, default=0.0, help="Degrees. Cut plane angle positive in ccw direction from X axis.")
parser.add_argument('-strain', action='store_true', help="Inputs are strain values.")
parser.add_argument('--E', type=float, default=0.0, help="Elastic Modulus.")
parser.add_argument('--v', type=float, default=0.0, help="Poisson's Ratio.")
parser.add_argument('-viz', action='store_true', help="Visualize Stress Element.")
args = parser.parse_args()

sigma_x = args.sx
sigma_y = args.sy
tau_xy = args.txy
theta_deg = args.deg  # degrees, ccw is positive

if args.strain:
    tau_xy = tau_xy/2  # MOHR CIRCLE USES HALF SHEAR FOR STRAIN

normal_symbol = "ε" if args.strain else "σ"
shear_symbol = "γ" if args.strain else "τ"
value_unit = "Strain" if args.strain else "Stress"

if args.strain:
    print("INPUTS ARE STRAIN")
else:
    print("INPUTS ARE STRESS")

print(f"{normal_symbol}_x = {sigma_x:e}")
print(f"{normal_symbol}_y = {sigma_y:e}")
print(f"{shear_symbol}_xy = {tau_xy:e}")
print(f"θ = {theta_deg}")

tau_max_in_plane = np.sqrt( (0.5*(sigma_x-sigma_y))**2 + tau_xy**2 )
sigma_avg = 0.5*(sigma_x+sigma_y)

if args.deg != 0.0:
    print(f"\n{value_unit} transformation")
    print(f"{value_unit} at {theta_deg} degrees.")
    sigma_x_prime = (0.5*(sigma_x+sigma_y)) + (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_deg))) + (tau_xy*np.sin(2*np.radians(theta_deg)))
    sigma_y_prime = (0.5*(sigma_x+sigma_y)) - (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_deg))) - (tau_xy*np.sin(2*np.radians(theta_deg)))
    tau_xy_prime = (-0.5*(sigma_x-sigma_y)*np.sin(2*np.radians(theta_deg)) + (tau_xy*np.cos(2*np.radians(theta_deg))))
    print(f"{normal_symbol}_x' = {sigma_x_prime:e}")
    print(f"{normal_symbol}_y' = {sigma_y_prime:e}")
    print(f"{shear_symbol}_x'y' = {tau_xy_prime:e}")


print(f"\nPrincipal {value_unit}, Plane with no Shear {value_unit}.")

# print("in_top :", (2*tau_xy))
# print("in_bot :", (sigma_x-sigma_y))
# print("in_atan :", (2*tau_xy) / (sigma_x-sigma_y)) 
theta_p1 = 0.5 * np.degrees( np.arctan2(  (2*tau_xy) , (sigma_x-sigma_y) ) )
# theta_p1 = 0.5 * np.degrees( np.arctan( (2*tau_xy) / (sigma_x-sigma_y) ) )

# tan(2a) = tan(2a+180) = b
# hence, a = arctan(a)/2 - 90
theta_p2 = ( theta_p1 - 90 )
while theta_p2 < -90:
    theta_p2 += 180

print(f"Orientation of Principal Planes of {value_unit}")

print(f"θ_p1={theta_p1:.6f} deg")
sigma_x_prime_p1 = (0.5*(sigma_x+sigma_y)) + (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_p1))) + (tau_xy*np.sin(2*np.radians(theta_p1)))
sigma_y_prime_p1 = (0.5*(sigma_x+sigma_y)) - (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_p1))) - (tau_xy*np.sin(2*np.radians(theta_p1)))
tau_xy_prime1 = (-0.5*(sigma_x-sigma_y)*np.sin(2*np.radians(theta_p1)) + (tau_xy*np.cos(2*np.radians(theta_p1))))
print(f"p1: {normal_symbol}_x'={sigma_x_prime_p1:e}, {normal_symbol}_y'={sigma_y_prime_p1:e}, {shear_symbol}_x'y'={tau_xy_prime1:e}")

print(f"θ_p2={theta_p2:.6f} deg")
sigma_x_prime_p2 = (0.5*(sigma_x+sigma_y)) + (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_p2))) + (tau_xy*np.sin(2*np.radians(theta_p2)))
sigma_y_prime_p2 = (0.5*(sigma_x+sigma_y)) - (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_p2))) - (tau_xy*np.sin(2*np.radians(theta_p2)))
tau_xy_prime2 = (-0.5*(sigma_x-sigma_y)*np.sin(2*np.radians(theta_p2)) + (tau_xy*np.cos(2*np.radians(theta_p2))))
print(f"p2: {normal_symbol}_x'={sigma_x_prime_p2:e}, {normal_symbol}_y'={sigma_y_prime_p2:e}, {shear_symbol}_x'y'={tau_xy_prime2:e}")


print(f"\nMaximum In-Plane Shear {value_unit}")

theta_s1 = 0.5 * np.degrees( np.arctan2( (sigma_x-sigma_y) , (-2*tau_xy) ) )
# theta_s1 = 0.5 * np.degrees( np.arctan( (sigma_x-sigma_y) / (-2*tau_xy) ) )

# tan(2a) = tan(2a+180) = b
# hence, a = arctan(a)/2 - 90
theta_s2 = ( theta_s1 - 90 )
while theta_s2 < -90:
    theta_s2 += 180

print(f"θ_s1={theta_s1:.6f} deg. θ_s2={theta_s2:.6f} deg.")

if theta_s1 > 0:
    tau_max_in_plane *= -1


print(f"{'γ/2' if args.strain else 'τ'} max in plane = {tau_max_in_plane:e}. **Sign matters.")
print("Use '-viz' argument to visualize the In-Plane Shear.")
print(f"{normal_symbol} avg = {sigma_avg:e}")


print(f"\nAbsolute Maximum Shear {value_unit}")

sigma_1 = sigma_avg + abs(tau_max_in_plane)
sigma_2 = sigma_avg - abs(tau_max_in_plane)

# if abs(sigma_1) < abs(sigma_2):
#     sigma_1, sigma_2 = sigma_2, sigma_1

print(f"{normal_symbol}_1={sigma_1:e}, {normal_symbol}_2={sigma_2:e}")

if sigma_1 * sigma_2 > 0 :  # same sign
    print(f"{normal_symbol}_1 and {normal_symbol}_2 have the same sign.")
    print(f"{normal_symbol}_3 is out of plane.")
    if abs(sigma_1) > abs(sigma_2):
        # print(f"{normal_symbol}_1 is larger.")
        tau_abs_max = 0.5*sigma_1
    else:
        # print(f"{normal_symbol}_2 is larger.")
        tau_abs_max = 0.5*sigma_2
else:  # opposite sign
    print(f"{normal_symbol}_1 and {normal_symbol}_2 have the opposite sign.")
    print(f"maximum shear equals {'γ/2' if args.strain else 'τ'} max in-plane.")
    tau_abs_max = 0.5*(sigma_1-sigma_2)

print(f"{'γ/2' if args.strain else 'τ'} abs max = {tau_abs_max:e}")


print("\nMohr's Circle")
print(f"X point = ({sigma_x:e}, {tau_xy:e})")
print("Rotate X to the closest X axis.")
print(f"2 θp = {np.degrees(np.arctan2(-tau_xy, sigma_x-sigma_avg))} from the σ+ axis.")
print(f"Y point = ({sigma_y:e}, {-tau_xy:e})")
print(f"Center ({normal_symbol} avg) = ({sigma_avg:e}, 0)")
print(f"Radius ({'γ/2' if args.strain else 'τ'} max) = {abs(tau_max_in_plane):e}")


if not args.strain:
    print("\n*** REMINDERS ***")
    print(" |θp| + |θs| = 45 degrees.")
    print(" σ3 is always at zero. Out of plane constant.")

if args.strain:
    if args.E and args.v:
        print("\nPrincipal Stress from Strain, Hooke's Law for plane stress")
        # calcuate the sigma_# from the stresses, bc rn sigma_# is the strain input
        ev = args.E/(1-(args.v*args.v))
        sigma_1, sigma_2 = ev*(sigma_1 + args.v*sigma_2), ev*(sigma_2 + args.v*sigma_1)
        print(f"σ_1={sigma_1:e}")
        print(f"σ_2={sigma_2:e}")

        if sigma_1 * sigma_2 > 0 :  # same sign
            print("σ_1 and σ_2 have the same sign.")
            print("σ_3 is out of plane.")
            if abs(sigma_1) > abs(sigma_2):
                # print("σ_1 is larger.")
                tau_abs_max = 0.5*sigma_1
            else:
                # print("σ_2 is larger.")
                tau_abs_max = 0.5*sigma_2
        else:  # opposite sign
            print("σ_1 and σ_2 have the opposite sign.")
            print("maximum shear equals τ max in-plane.")
            tau_abs_max = 0.5*(sigma_1-sigma_2)
        print(f"τ abs max = {tau_abs_max:e}")


if args.viz:
    from PIL import Image, ImageDraw, ImageFont

    # Size and position
    size = (1200, 1000)
    scale = 90
    font_size = int(scale / 4.8)
    # print("font_size :", font_size)
    font_path = "fonts/InputSans-Regular.ttf"


    step_scale = 4.4
    step_x, step_y = size[0]/step_scale, size[1]/step_scale

    input_center = (step_x, 2*step_y)
    principal_center = (3*step_x, step_y)
    shear_center = (3*step_x, 3*step_y)

    # Colors
    bkg_color = (192, 192, 192)
    arrow_color = (222,184,135)
    arrow_outline = (0, 0, 0)
    square_color = (128, 128, 128)
    square_outline = (0, 0, 0)
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
        draw.polygon(xy=tuple([tuple(x) for x in input_square]), fill=square_color, outline=square_outline)

        text_move = 1.75

        # shear arrows
        if txy != 0:
            arr = rotate_2d(arrow_points(scale), 90+angle) + rotate_2d([scale*1.25, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            arr = rotate_2d(arrow_points(scale), -90+angle) + rotate_2d([-scale*1.25, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            arr = rotate_2d(arrow_points(scale), 0+angle) + rotate_2d([scale*1.25, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            arr = rotate_2d(arrow_points(scale), 180+angle) + rotate_2d([-scale*1.25, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)

            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], 225+angle, center)), text="txy={:e}".format(txy), fill=text_color, anchor="mm", font=font)

        # normal arrows
        if sx != 0:
            arr = rotate_2d(arrow_points(scale*0.6), 0+angle) + rotate_2d([scale*1.6, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            arr = rotate_2d(arrow_points(scale*0.6), 180+angle) + rotate_2d([-scale*1.6, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], 180+angle, center)), text="sx={:e}".format(sx), fill=text_color, anchor="mm", font=font)

        if sy != 0:
            arr = rotate_2d(arrow_points(scale*0.6), 90+angle) + rotate_2d([scale*1.6, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            arr = rotate_2d(arrow_points(scale*0.6), -90+angle) + rotate_2d([-scale*1.6, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], -90+angle, center)), text="sy={:e}".format(sy), fill=text_color, anchor="mm", font=font)

        # Title
        draw.text(xy=tuple(rotate_2d([0, -scale*2.35], 0, center)), text="{}, θ={:.2f}.".format(title, angle), fill=text_color, anchor="mm", font=font)

    # Draw element
    draw_element(center=input_center, scale=scale, title="Input", angle=args.deg, sx=sigma_x, sy=sigma_y, txy=tau_xy)

    # Draw pricinple axis, no shear
    draw_element(center=principal_center, scale=scale, title="Principal Plane (No Shear)", angle=theta_p1, sx=sigma_x_prime_p1, sy=sigma_y_prime_p1)

    # Draw shear axis
    draw_element(center=shear_center, scale=scale, title="Shear Plane (Maximum Shear)", angle=theta_s1, sx=sigma_avg, sy=sigma_avg, txy=tau_max_in_plane)


    im.show()


