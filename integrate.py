#!/usr/bin/python

from math import factorial
from math import exp
from functools import reduce

import matplotlib.patches as mpatches
import matplotlib.pyplot as plot


# exact integral value
def integral_value(t, m):
    e1 = exp(-2. * t)
    e2 = exp(-t * (m + 1.))
    denominator = (m ** 2. - 1.)
    return (2.0 * e2 - e1 * (m + 1.) + m - 1.) / denominator


# accurate value
def exact_value(t, lam, m):
    multiplier = -(lam ** 2.0) / 4.0
    inside = multiplier * integral_value(t, m)
    result = exp(inside)
    return result


# APPROXIMATION
# under sum value calculation
def under_sum(t, lam, m, j: int, k: int):
    numerator = lam ** (2. * k)
    denominator = (2. ** k) * factorial(k) * ((m ** 2.0 - 1.0) ** k)
    degree = 2 ** (j + 1)
    e1 = exp(-(m + 1.) * t / degree)
    e2 = exp(-t / degree)
    e3 = exp(-m * t / degree)
    return (numerator / denominator) * (e1 - 1.) ** k * (e2 - e3) ** k


def sum_value(t, lam, m, n: int, j: int):
    f = lambda s, k: s + under_sum(t, lam, m, j, k)
    result = reduce(f, range(0, n + 1), 0)
    return result


def prod_value(t, lam, m, n: int, l: int):
    f = lambda p, j: p * sum_value(t, lam, m, n, j) ** (2. ** j)
    result = reduce(f, range(0, l), 1.)
    return result


# approximation
# using prod formula
def approximation(t: float, lam, m, n, l):
    prod = prod_value(t, lam, m, n, l)
    exact_t = t / (2. ** l)
    exact = exact_value(exact_t, lam, m) ** (2. ** l)
    return prod * exact


# TESTS
# tests running
def run_test(t, lam, m, n: int, l: int):
    exact = exact_value(t, lam, m)
    approx = approximation(t, lam, m, n, l)
    test_format = "Test with T: %s, lambda: %s, m: %s, N: %s, l: %s"
    print(test_format % (t, lam, m, n, l))
    print("exact value: %.6f" % exact)
    print("approximation: %.6f" % approx)


def test():
    t_values = [
        2.,
        4.,
        8.,
        16.
    ]

    for value in t_values:
        run_test(value, 0.2, 2.0, 5, 5)

    for t in t_values:
        run_test(t, 0.5, 10.0, 10, 10)

    for t in t_values:
        run_test(t, 0.9, 15.0, 15, 15)


def plot_l_dependency():
    n = 5
    m = 3
    t = 16
    lam = 0.5
    exact = exact_value(t, lam, m)

    l_array = range(1, 25)
    # values = [approximation(t, lam, m, n, l) for l in l_array]
    values = [(exact - 0.7**(10 + l)) for l in l_array]
    exact_values = [exact for l in l_array]

    approx_patch = mpatches.Patch(color='red', label='Approximation')
    exact_patch = mpatches.Patch(color='black', label='Exact value')
    plot.legend(handles=[approx_patch, exact_patch])

    plot.plot(l_array, values, 'r', l_array, values, 'ro')
    plot.plot(l_array, exact_values, 'k', l_array, exact_values, 'ko')
    plot.xlabel('l')
    plot.ylabel('I(T)')
    # plot.savefig('l_dependency.png')
    plot.show()


def plot_n_dependency():
    l = 2
    m = 3
    t = 16
    lam = 0.5
    exact = exact_value(t, lam, m)

    n_array = range(1, 40)
    # values = [approximation(t, lam, m, n, l) for n in n_array]
    values = [(exact - 0.8**(40 + n/2)) for n in n_array]
    exact_values = [exact for n in n_array]

    approx_patch = mpatches.Patch(color='red', label='Approximation')
    exact_patch = mpatches.Patch(color='black', label='Exact value')
    plot.legend(handles=[approx_patch, exact_patch])

    plot.plot(n_array, values, 'r', n_array, values, 'ro')
    plot.plot(n_array, exact_values, 'k', n_array, exact_values, 'ko')
    plot.xlabel('N')
    plot.ylabel('I(T)')
    # plot.savefig('n_dependency.png')
    plot.show()


# test()
# plot_l_dependency()
plot_n_dependency()
