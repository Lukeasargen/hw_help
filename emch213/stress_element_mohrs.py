import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Shear Tranforms and Mohr's Circle")
parser.add_argument('--sx', type=float, default=0.0, metavar="Stress in X axis")
parser.add_argument('--sy', type=float, default=0.0, metavar="Stress in Y axis")
parser.add_argument('--txy', type=float, default=0.0, metavar="Shear Stress.")
parser.add_argument('--deg', type=float, default=0.0, metavar="Degrees. Cut plane angle positive in ccw direction from X axis.")
parser.add_argument('--v', action='store_true')
args = parser.parse_args()

sigma_x = args.sx
sigma_y = args.sy
tau_xy = args.txy
theta_deg = args.deg  # degrees, ccw is positive
visualize = args.v

print("σ x = {}".format(sigma_x))
print("σ x = {}".format(sigma_y))
print("τ xy = {}".format(tau_xy))
print("θ = {}".format(theta_deg))


print("\nStress transformation")

tau_max_in_plane = np.sqrt( (0.5*(sigma_x-sigma_y))**2 + tau_xy**2 )
sigma_avg = 0.5*(sigma_x+sigma_y)


print("Stresses at {} degrees.".format(theta_deg))
sigma_x_prime = (0.5*(sigma_x+sigma_y)) + (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_deg))) + (tau_xy*np.sin(2*np.radians(theta_deg)))
sigma_y_prime = (0.5*(sigma_x+sigma_y)) - (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_deg))) - (tau_xy*np.sin(2*np.radians(theta_deg)))
tau_xy_prime = (-0.5*(sigma_x-sigma_y)*np.sin(2*np.radians(theta_deg)) + (tau_xy*np.cos(2*np.radians(theta_deg))))
print("σx' = {}".format(sigma_x_prime))
print("σy' = {}".format(sigma_y_prime))
print("τx'y' = {}".format(tau_xy_prime))


print("\nPrincipal Stress")
print("Plane with no shear stress.")

# print("in_top :", (2*tau_xy))
# print("in_bot :", (sigma_x-sigma_y))
# print("in_atan :", (2*tau_xy) / (sigma_x-sigma_y)) 
theta_p1_temp = 0.5 * np.degrees( np.arctan2(  (2*tau_xy) , (sigma_x-sigma_y) ) )
# theta_p1_temp = 0.5 * np.degrees( np.arctan( (2*tau_xy) / (sigma_x-sigma_y) ) )

# tan(2a) = tan(2a+180) = b
# hence, a = arctan(a)/2 - 90
theta_p2_temp = ( theta_p1_temp - 90 )
while theta_p2_temp < -90:
    theta_p2_temp += 180

if abs(theta_p1_temp) > abs(theta_p2_temp):
    theta_p1 = theta_p2_temp
    theta_p2 = theta_p1_temp
else:
    theta_p1 = theta_p1_temp
    theta_p2 = theta_p2_temp

print("Orientation of Principal Planes of Stress")

print("θ p1={:.6f}".format(theta_p1))
sigma_x_prime_p1 = (0.5*(sigma_x+sigma_y)) + (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_p1))) + (tau_xy*np.sin(2*np.radians(theta_p1)))
sigma_y_prime_p1 = (0.5*(sigma_x+sigma_y)) - (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_p1))) - (tau_xy*np.sin(2*np.radians(theta_p1)))
tau_xy_prime1 = (-0.5*(sigma_x-sigma_y)*np.sin(2*np.radians(theta_p1)) + (tau_xy*np.cos(2*np.radians(theta_p1))))
print("σx' p1={:.6f}, σy' p1={:.6f}, τx'y' p1={:.6f}".format(sigma_x_prime_p1, sigma_y_prime_p1, tau_xy_prime1))

print("θ p2={:.6f}".format(theta_p2))
sigma_x_prime_p2 = (0.5*(sigma_x+sigma_y)) + (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_p2))) + (tau_xy*np.sin(2*np.radians(theta_p2)))
sigma_y_prime_p2 = (0.5*(sigma_x+sigma_y)) - (0.5*(sigma_x-sigma_y)*np.cos(2*np.radians(theta_p2))) - (tau_xy*np.sin(2*np.radians(theta_p2)))
tau_xy_prime2 = (-0.5*(sigma_x-sigma_y)*np.sin(2*np.radians(theta_p2)) + (tau_xy*np.cos(2*np.radians(theta_p2))))
print("σx' p2={:.6f}, σy' p2={:.6f}, τx'y' p2={:.6f}".format(sigma_x_prime_p2, sigma_y_prime_p2, tau_xy_prime2))


print("\nMaximum In-Plane Shear Stress")

theta_s1_temp = 0.5 * np.degrees( np.arctan2( (sigma_x-sigma_y) , (-2*tau_xy) ) )
# theta_s1_temp = 0.5 * np.degrees( np.arctan( (sigma_x-sigma_y) / (-2*tau_xy) ) )

# tan(2a) = tan(2a+180) = b
# hence, a = arctan(a)/2 - 90
theta_s2_temp = ( theta_s1_temp - 90 )
while theta_s2_temp < -90:
    theta_s2_temp += 180

if abs(theta_s1_temp) > abs(theta_s2_temp):
    theta_s1 = theta_s2_temp
    theta_s2 = theta_s1_temp
else:
    theta_s1 = theta_s1_temp
    theta_s2 = theta_s2_temp

print("θ s1={:.6f}, θ s2={:.6f}".format(theta_s1, theta_s2))

if theta_s1 > 0:
    tau_max_in_plane *= -1

print("τ max in plane = {}. **Sign matters.\nUse '--v' argument to visualize the In-Plane Shear.".format(tau_max_in_plane))
print("σ avg = {}".format(sigma_avg))


print("\nAbsolute Maximum Shear Stress")

sigma_1_temp = sigma_avg + abs(tau_max_in_plane)
sigma_2_temp = sigma_avg - abs(tau_max_in_plane)

if abs(sigma_1_temp) > abs(sigma_2_temp):
    sigma_1 = sigma_1_temp
    sigma_2 = sigma_2_temp
else:
    sigma_1 = sigma_2_temp
    sigma_2 = sigma_1_temp

print("σ1={:.6f}, σ2={:.6f}".format(sigma_1, sigma_2))

if sigma_1 * sigma_2 > 0 :  # same sign
    print("σ1 and σ2 have the same sign.")
    if abs(sigma_1) > abs(sigma_2):
        tau_abs_max = 0.5*sigma_1
    else:
        tau_abs_max = 0.5*sigma_2
else:  # opposite sign
    print("σ1 and σ2 have the opposite sign.")
    tau_abs_max = 0.5*(sigma_1-sigma_2)

print("τ abs max = {}".format(tau_abs_max))


print("\nMohr's Circle")


print("X point = ({:.6f}, {:.6f})".format(sigma_x, -tau_xy))
print("Rotate X to the closest X axis.")
print("2 θp = {}".format( np.degrees( np.arctan2( -tau_xy, sigma_x-sigma_avg )) ) )
print("Y point = ({:.6f}, {:.6f})".format(sigma_y, tau_xy))
print("Center (σ avg) = ({:.6f}, 0)".format(sigma_avg))
print("Radius (τ max) = {:.6f}".format(tau_max_in_plane))


print("\n*** REMINDERS ***")
print(" |θp| + |θs| = 45 degrees.")
print(" σ3 is always at zero.")


if visualize:
    from PIL import Image, ImageDraw, ImageFont

    # Size and position
    size = (1000, 1000)
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

            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], 225+angle, center)), text="txy={:.2f}".format(txy), fill=text_color, anchor="mm", font=font)

        # normal arrows
        if sx != 0:
            arr = rotate_2d(arrow_points(scale*0.6), 0+angle) + rotate_2d([scale*1.6, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            arr = rotate_2d(arrow_points(scale*0.6), 180+angle) + rotate_2d([-scale*1.6, 0], angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], 180+angle, center)), text="sx={:.2f}".format(sx), fill=text_color, anchor="mm", font=font)

        if sy != 0:
            arr = rotate_2d(arrow_points(scale*0.6), 90+angle) + rotate_2d([scale*1.6, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            arr = rotate_2d(arrow_points(scale*0.6), -90+angle) + rotate_2d([-scale*1.6, 0], 90+angle, center)
            draw.polygon(xy=tuple([tuple(x) for x in arr]), fill=arrow_color, outline=arrow_outline)
            draw.text(xy=tuple(rotate_2d([-scale*text_move, 0], -90+angle, center)), text="sy={:.2f}".format(sy), fill=text_color, anchor="mm", font=font)

        # Title
        draw.text(xy=tuple(rotate_2d([-scale*2.25, -scale*2.25], 0, center)), text="{}, θ={:.2f}.".format(title, angle), fill=text_color, anchor="mm", font=font)

    # Draw element
    draw_element(center=input_center, scale=scale, title="Input", sx=sigma_x, sy=sigma_y, txy=tau_xy)

    # Draw pricinple axis, no shear
    draw_element(center=principal_center, scale=scale, title="Principal Plane (No Shear)", angle=theta_p1, sx=sigma_x_prime_p1, sy=sigma_y_prime_p1)

    # Draw shear axis
    draw_element(center=shear_center, scale=scale, title="Shear Plane (Maximum Shear)", angle=theta_s1, sx=sigma_avg, sy=sigma_avg, txy=tau_max_in_plane)


    im.show()


