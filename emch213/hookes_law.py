import argparse

parser = argparse.ArgumentParser(description="Shear Tranforms and Mohr's Circle")
parser.add_argument('--v', type=float, default=0.0, required=True, metavar="Poisson Ratio")
parser.add_argument('--e', type=float, default=0.0, required=True, metavar="Elastic Modulus")
parser.add_argument('--sx', type=float, default=0.0, metavar="Stress Z Axis")
parser.add_argument('--sy', type=float, default=0.0, metavar="Stress Z Axis")
parser.add_argument('--sz', type=float, default=0.0, metavar="Stress Z Axis")
args = parser.parse_args()


sigma_x = args.sx
sigma_y = args.sy
sigma_z = args.sz
poisson_ratio = args.v
elastic_modulus = args.e

print("σx = {:.6}".format(sigma_x))
print("σy = {:.6}".format(sigma_y))
print("σz = {:.6}".format(sigma_z))
print("v (Poisson Ratio) = {}".format(poisson_ratio))
print("E (Elastic Modulus) = {:.6}".format(elastic_modulus))


strain_x = (1/elastic_modulus) * ( sigma_x - ( poisson_ratio * (sigma_y + sigma_z) ) )
strain_y = (1/elastic_modulus) * ( sigma_y - ( poisson_ratio * (sigma_z + sigma_x) ) )
strain_z = (1/elastic_modulus) * ( sigma_z - ( poisson_ratio * (sigma_x + sigma_y) ) )
print("Strain X = {:.6}".format(strain_x))
print("Strain Y = {:.6}".format(strain_y))
print("Strain Z = {:.6}".format(strain_z))


shear_modulus = elastic_modulus / (2 * (1+poisson_ratio) )
print("G (Shear Modulus Calculated) = {:.4}".format(shear_modulus))
