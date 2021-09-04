import numpy as np
from scipy.optimize import minimize


esp = 0.002e-3
D = 0.05
Re = 201050

def colebrook(f):
    A = (esp/D)/3.2
    B = 2.51/(Re*np.sqrt(f))
    c = -2.0*np.log10(A+B) - (1/np.sqrt(f))
    print(f, c)
    return c

def objective(f):
    return colebrook(f)

def eq(f):
    return colebrook(f)

sol = minimize(objective, [1e-4], constraints={'type': 'eq', 'fun': eq})

print(sol)
