import numpy as np
import matplotlib.pyplot as plt
'from math import cos, sin, pi'

def f1(x):
    return np.sin(np.cos(x))

def f2(x):
    return np.abs(np.abs(x) - 1)

def equal_nodes(a, b, n):
    return np.linspace(a, b, n)

def chebyshev_nodes(a, b, n):
    k = np.arange(1, n + 1)
    return (a + b) / 2 + (b - a) / 2 * np.cos((2 * k - 1) * np.pi / (2 * n))

def compute_newton_coefficients(x, f):
    n = len(x)
    F = np.zeros((n, n))
    F[:,0] = f(x)

    for j in range(1, n):
        for i in range(n - j):
            F[i, j] = (F[i + 1, j - 1] - F[i, j - 1]) / (x[i + j] - x[i])

    return F[0]

def evaluate_newton_polynomial(x, nodes, coefficients):
    n = len(nodes)
    p = coefficients[-1]
    for i in range(n - 2, -1, -1):
        p = p * (x - nodes[i]) + coefficients[i]
    return p

def max_interpolation_error(f, nodes, coefficients, a, b, num_points=1000):
    x_values = np.linspace(a, b, num_points)
    y_values = f(x_values)
    interpolated_values = evaluate_newton_polynomial(x_values, nodes, coefficients)
    errors = np.abs(interpolated_values - y_values)
    return np.max(errors)

# Функция для построения графика интерполяции
def plot_interpolation(f, nodes, coefficients, a, b, label):
    x_values = np.linspace(a, b, 1000)
    y_values = f(x_values)
    interpolated_values = evaluate_newton_polynomial(x_values, nodes, coefficients)

    plt.plot(x_values, y_values, label="Original Function")
    plt.plot(x_values, interpolated_values, label=label)
    plt.scatter(nodes, f(nodes), color='red', label='Interpolation nodes')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)

def plot_interpolation2(f, nodes, coefficients, a, b, label):
    x_values = np.linspace(a, b, 1000)
    y_values = f(x_values)
    interpolated_values = evaluate_newton_polynomial(x_values, nodes, coefficients)
    plt.plot(x_values, interpolated_values, label=label)
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)

# Определение интервала
a, b = -3, 3

# Задание количества узлов интерполяции для каждой функции
n_values = [3, 10, 20]

# Создание общего графического окна для всех графиков
plt.figure(figsize=(12, 8))


# Интерполяция для обеих функций с использованием равномерных и Чебышёвских узлов
for i, (f, label) in enumerate(zip([f1, f2], ["f1(x) = sin(cos(x))", "f2(x) = | |x| - 1 | |"])):
    for n in n_values:
        plt.subplot(2, len(n_values), i * len(n_values) + n_values.index(n) + 1)

        # Равномерно распределенные узлы
        equal_x = equal_nodes(a, b, n)
        equal_coefficients = compute_newton_coefficients(equal_x, f)
        plot_interpolation(f, equal_x, equal_coefficients, a, b, "Equal nodes")

        # Узлы Чебышёва
        chebyshev_x = chebyshev_nodes(a, b, n)
        chebyshev_coefficients = compute_newton_coefficients(chebyshev_x, f)
        plot_interpolation2(f, chebyshev_x, chebyshev_coefficients, a, b, "Chebyshev nodes")

        plt.title(f'{label}, n={n}')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid(True)

plt.tight_layout()
plt.show()

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
        plot_interpolation(f, equal_x, equal_coefficients, a, b, "Равномерные узлы")

        # Узлы Чебышёва
        chebyshev_x = chebyshev_nodes(a, b, n)
        chebyshev_coefficients = compute_newton_coefficients(chebyshev_x, f)
        chebyshev_error = max_interpolation_error(f, chebyshev_x, chebyshev_coefficients, a, b)
        print("Погрешность интерполяции для узлов Чебышёва:", chebyshev_error)
        plot_interpolation(f, chebyshev_x, chebyshev_coefficients, a, b, "Узлы Чебышёва")

        print("=======================")


