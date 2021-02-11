import numpy as np
from scipy.optimize import minimize

def objective(X):
    x, y, z = X
    return x**2 + y**2 + z**2

def eq(X):
    x, y, z = X
    return x + 6*y + 7*z - 6

sol = minimize(objective, [1, -0.5, 0.5], constraints={'type': 'eq', 'fun': eq})

print(sol)

