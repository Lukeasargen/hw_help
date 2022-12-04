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

gamma = 1.4
M1 = np.array([2.38])  # mach number
theta = 11.9  # degrees
beta_weak = beta_from_M_theta(beta_guess=10, M=M1, theta=theta, gamma=gamma)[0]
print(f"{beta_weak = }")
beta_strong = beta_from_M_theta(beta_guess=80, M=M1, theta=theta, gamma=gamma)[0]
print(f"{beta_strong = }")
