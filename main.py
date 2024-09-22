import sympy as sp
from sympy import *
from sympy.utilities.lambdify import lambdify
x = symbols('x')

def Lagrange_i(data_points, i):
    """
    :param data_points: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th lagrange function
    """
    li, xi = 1, data_points[i][0]

    for j in range(len(data_points)):
        xj = data_points[j][0]

        if j != i:
            li *= (x - xj) / (xi - xj)
    return li


def Lagrange_squared(data_points, i):
    """
    :param data_points: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th squared lagrange non-function
    """
    l2 = Lagrange_i(data_points, i) ** 2
    return l2

def reduce_func(func):
    reduced_func = func
    for val in preorder_traversal(func):
        if isinstance(val, Float):
            reduced_func = reduced_func.subs(val, round(val, 1))
    return reduced_func

def Lagrange_Derivative(data_points, i):
    """
    :param data_points: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th Derivative lagrange function
    """
    dl = Lagrange_i(data_points, i).diff(x)
    return lambdify(x, dl)


def Hermite_Interpolation(data_points):
    """
    :param data_points: list of trios (x, f(x), f'(x))
    :return: the hermite function
    """
    print("\nfor the Data points:\n" + str(data_points) + "\n")
    hermite_function = 0
    for i in range(len(data_points)):
        xi, yi, mi = data_points[i]
        la_squared = Lagrange_squared(data_points, i)
        la_derivative = Lagrange_Derivative(data_points, i)
        hermite_function += (1 - 2*(x - xi) * la_derivative(xi)) * la_squared * yi + (x - xi)*la_squared * mi
    h = hermite_function
    dh = hermite_function.diff(x)
    reduced_h = reduce_func(simplify(h))
    reduced_dh = reduce_func(simplify(dh))
    h = lambdify(x, h)
    dh = lambdify(x, dh)
    for xi, yi, mi in data_points:
        print("x=" + str(xi) + " f(x)=" + str(yi) + " f'(x)=" + str(mi) + " H(x)=" + str(round(h(xi), 2)) + " H'(x)=" + str(round(dh(xi), 2)))
    return reduced_h, reduced_dh


def Drive():
    """
    getting user data and calculate H(x)
    :return: None
    """
    print("please put you Data points as follows:")
    data_points = list()
    while True:
        xi = float(input("xi: "))
        yi = float(input("yi: "))
        mi = float(input("mi: "))
        data_points.append((xi, yi, mi))

        stop = input("Continue? N - no, any other key - Yes ")
        if stop == "n" or stop == "N":
            break

    h, dh = Hermite_Interpolation(data_points)
    print("H(x) = " + str(h))
    print("H'(x) = " + str(dh))
    h = lambdify(x, h)
    dh = lambdify(x, dh)
    input_var = int(input("\nEnter x value to calculate H(x) and H'(x):\n"))
    print(f'H({input_var}) = {round(h(input_var), 2)}' + ", H'" + f'({input_var}) = {round(dh(input_var), 2)}')


Drive()

















