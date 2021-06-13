import math

import sympy as sp
from sympy.utilities.lambdify import lambdify

def Newton_Raphson(pol, little_range, epsilon):
    f = lambdify(x, pol)
    df = lambdify(x, sp.diff(pol, x))
    x1 = sum(little_range) / 2
    x2 = x1 - f(x1) / df(x1)
    counter = 1

    while abs(x2 - x1) > epsilon:
        x1 = x2
        x2 = (x1 - f(x1) / df(x1))
        counter += 1

    return x2, counter


def secant_method(pol, little_range, epsilon):
    f = lambdify(x, pol)
    x1, x2 = little_range
    x3 = (x1*f(x2) - x2*f(x1)) / (f(x2) - f(x1))
    counter = 1

    while abs(x3 - x2) > epsilon:
        counter += 1
        x1 = x2
        x2 = x3
        x3 = (x1*f(x2) - x2*f(x1)) / (f(x2) - f(x1))

    return x3, counter

x = sp.symbols('x')
user_func = (math.sin(2*x**3 + 5*x**2 - 6))/2*math.e**(-2*x)
user_range = (-1, 1.5)
user_epsilon = 0.0001
user_step = 0.1

def Driver(f)