import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve


def beta_from_M_theta(beta_guess, M, theta, gamma=1.4):
    def f(beta_guess, M, theta, gamma):
        b_rad = np.radians(beta_guess)
        top = (M*np.sin(b_rad))**2 - 1
        bot = M*M*(gamma + np.cos(2*b_rad)) + 2
        return np.tan(np.radians(theta)) - 2*(top/bot)/np.tan(b_rad)
    beta_guess = np.array([beta_guess]*len(M))
    beta = fsolve(f, beta_guess, args=(M, theta, gamma))
    return beta

def PM_function(M, gamma=1.4):
    part1 = np.sqrt( (gamma+1)/(gamma-1) )
    part2 = np.arctan( np.sqrt( (M*M - 1)*(gamma-1)/(gamma+1) ) )
    part3 = np.arctan( np.sqrt( M*M - 1) )
    return np.degrees(part1*part2 - part3)

def PM_function_from_nu(M_guess, nu, gamma=1.4):
    def f(M_guess, nu, gamma):
        return nu - PM_function(M_guess, gamma)
    M_guess = np.array([M_guess]*len(nu))
    M = fsolve(f, M_guess, args=(nu, gamma))
    return M

gamma = 1.4
M1 = np.array([2.38])  # mach number

nu = 33.8
M_out = PM_function_from_nu(M_guess=1.5, nu=[nu], gamma=gamma)[0]
print(f"{M_out = }")

