
from sympy import sin
from sympy import exp
import sympy as sp
from sympy.utilities.lambdify import lambdify

def Drive(pol, big_range, epsilon, step):

    print("f(x) = ", pol, ",Range = ", big_range, ", epsilon = ", epsilon)
    print("---------------Secant Method---------------")
    Solver(pol, big_range, epsilon, step, secant_method)
    print("---------------Newton_Raphson---------------")
    Solver(pol, big_range, epsilon, step, Newton_Raphson)


def Solver(pol, big_range, epsilon, step, method):
    function_Solver(pol, big_range, epsilon, step, method)
    derivative_solver(pol, big_range, epsilon, step, method)

def function_Solver(pol, big_range, epsilon, step, method):
    f = lambdify(x, pol)
    left_bound, right_bound = big_range
    a, b = left_bound, left_bound + step

    while b <= right_bound:

        if f(a) * f(b) < 0:
            solution = method(pol, (a, b), epsilon)
            sol, iterations = solution
            if solution is not None:
                if abs(sol) < epsilon:
                    sol = 0
                print("x = " + str(sol) + ", number of iteration: " + str(iterations))
        a += step
        b += step


def derivative_solver(pol, big_range, epsilon, step, method):
    f = lambdify(x, pol)
    df = lambdify(x, sp.diff(pol, x))
    left_bound, right_bound = big_range
    a, b = left_bound, left_bound + step

    while b <= right_bound:

        if df(a) * df(b) < 0:

            if abs(df(b)) < epsilon or abs(df(a) < epsilon):
                solution = method(sp.diff(pol, x), (a - step, b + step), epsilon)
                sol, iterations = solution

            else:
                solution = method(sp.diff(pol, x), (a, b), epsilon)
                sol, iterations = solution

            if solution is not None and abs(f(sol)) < epsilon:
                if abs(sol) < epsilon:
                    sol = 0
                print("x = " + str(sol) + ", number of iterations : " + str(iterations))
                solution = None
        a += step
        b += step


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
user_func = (sin(2*x**3 + 5*x**2 - 6))/(2*exp(-2*x))
user_range = (-1, 1.5)
user_epsilon = 0.0001
user_step = 0.05

Drive(user_func, user_range, user_epsilon, user_step)