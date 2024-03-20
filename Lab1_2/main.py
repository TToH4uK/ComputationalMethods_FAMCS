import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi


# Определение функций f1(x) и f2(x)
def f1(x):
    return np.sin(np.cos(x))


def f2(x):
    return abs(abs(x) - 1)


# Генерация равномерно распределенных узлов на отрезке [a, b]
def equal_nodes(a, b, n):
    return np.linspace(a, b, n)


# Генерация узлов Чебышёва на отрезке [a, b]
def chebyshev_nodes(a, b, n):
    k = np.arange(1, n + 1)
    return (a + b) / 2 + (b - a) / 2 * np.cos((2 * k - 1) * np.pi / (2 * n))


# Вычисление коэффициентов интерполяционного многочлена Ньютона
def compute_newton_coefficients(x, f):
    n = len(x)
    F = np.zeros((n, n))
    F[:, 0] = f(x)

    for j in range(1, n):
        for i in range(n - j):
            F[i, j] = (F[i + 1, j - 1] - F[i, j - 1]) / (x[i + j] - x[i])

    return F[0]

# Вычисление значения интерполяционного многочлена Ньютона в точке x
def evaluate_newton_polynomial(x, nodes, coefficients):
    n = len(nodes)
    p = coefficients[-1]
    for i in range(n - 2, -1, -1):
        p = p * (x - nodes[i]) + coefficients[i]
    return p


# Вычисление максимальной погрешности интерполяции
def max_interpolation_error(f, nodes, coefficients, a, b, num_points=1000):
    x_values = np.linspace(a, b, num_points)
    y_values = f(x_values)
    interpolated_values = evaluate_newton_polynomial(x_values, nodes, coefficients)
    errors = np.abs(interpolated_values - y_values)
    return np.max(errors)


# Построение графика интерполяции с узлами равномерно распределенными и с узлами Чебышёва
def plot_interpolation_both(f, equal_nodes, equal_coefficients, chebyshev_nodes, chebyshev_coefficients, a, b, label):
    x_values = np.linspace(a, b, 1000)
    y_values = f(x_values)
    interpolated_values_equal = evaluate_newton_polynomial(x_values, equal_nodes, equal_coefficients)
    interpolated_values_chebyshev = evaluate_newton_polynomial(x_values, chebyshev_nodes, chebyshev_coefficients)

    plt.plot(x_values, y_values, label="Исходная функция")
    plt.plot(x_values, interpolated_values_equal, label="Равномерные узлы")
    plt.plot(x_values, interpolated_values_chebyshev, label="Узлы Чебышёва")
    plt.scatter(equal_nodes, f(equal_nodes), color='red', label='Узлы интерполяции (равномерные)')
    plt.scatter(chebyshev_nodes, f(chebyshev_nodes), color='green', label='Узлы интерполяции (Чебышёв)')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Интерполяция')
    plt.grid(True)
    plt.show()


# Определение интервала
a, b = -3, 3

# Задание количества узлов интерполяции для каждой функции
n_values = [3, 10, 20]

# Интерполяция для обеих функций с использованием равномерных и Чебышёвских узлов
for f, label in zip([f1, f2], ["f1(x) = sin(cos(x))", "f2(x) = | |x| - 1 | |"]):
    print("Функция:", label)
    print("=======================")

    for n in n_values:
        print(f"Интерполяция с {n} узлами")

        # Равномерно распределенные узлы
        equal_x = equal_nodes(a, b, n)
        equal_coefficients = compute_newton_coefficients(equal_x, f)
        equal_error = max_interpolation_error(f, equal_x, equal_coefficients, a, b)
        print("Погрешность интерполяции для равномерных узлов:", equal_error)

        # Узлы Чебышёва
        chebyshev_x = chebyshev_nodes(a, b, n)
        chebyshev_coefficients = compute_newton_coefficients(chebyshev_x, f)
        chebyshev_error = max_interpolation_error(f, chebyshev_x, chebyshev_coefficients, a, b)
        print("Погрешность интерполяции для узлов Чебышёва:", chebyshev_error)

        plot_interpolation_both(f, equal_x, equal_coefficients, chebyshev_x, chebyshev_coefficients, a, b, "Интерполяция")
        print("=======================")

