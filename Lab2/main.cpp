#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

const double eps = 1e-6; // Порог точности
const int kmax = 1000;   // Максимальное количество итераций
const int N = 5;

// Функция для генерации матрицы A
vector<vector<double>> GenerateMatrixA(int n) {
    vector<vector<double>> A(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        A[i][i] = 10.0 * sqrt(pow(i + 1, 3));
        for (int j = 0; j < n; j++) {
            if (i != j) {
                A[i][j] = 0.001 / pow(j + 1, 0.3);
            }
        }
    }

    return A;
}

// Функция для генерации вектора b
vector<double> GenerateVectorB(const vector<vector<double>>& A, int n) {
    vector<double> x_star(n);
    for (int i = 0; i < n; i++) {
        x_star[i] = i + 3; // Создаем вектор x* = (3, 4, 5, ...)
    }

    vector<double> b(n, 0.0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * x_star[j];
        }
    }

    return b;
}

bool IsStrictlyDiagonallyDominant(const vector<vector<double>>& A) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        double diagonal = abs(A[i][i]);
        double sum = 0.0;

        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += abs(A[i][j]);
            }
        }

        if (diagonal <= sum) {
            return false; // Не строго диагонально доминирующая матрица
        }
    }

    return true; // Строго диагонально доминирующая матрица
}

// Функция для вычисления суммы элементов A * X, исключая элемент с индексом excludeIndex
double SumAx(const vector<vector<double>>& A, const vector<double>& X, int excludeIndex) {
    double sum = 0.0;
    for (int i = 0; i < A.size(); i++) {
        if (i != excludeIndex) {
            sum += A[excludeIndex][i] * X[i];
        }
    }
    return sum;
}

// Метод Якоби
vector<double> JacobiMethod(const vector<vector<double>>& A, const vector<double>& b, int n) {
    vector<double> X(n, 0.0); // Начальное приближение

    int k = 0;
    while (k < kmax) {
        vector<double> X_new(n);
        for (int i = 0; i < n; i++) {
            X_new[i] = (b[i] - SumAx(A, X, i)) / A[i][i];
        }

        double maxDiff = 0.0;
        for (int i = 0; i < n; i++) {
            maxDiff = max(maxDiff, abs(X_new[i] - X[i]));
        }

        X = X_new;
        k++;

        if (maxDiff < eps) {
            cout << "Jacobi method: the required accuracy per iteration has been achieved on iter " << k << endl;
            return X;
        }
    }

    cout << "Jacobi method: the maximum number of iterations reached (" << kmax << ")" << endl;
    return X;
}

// Метод Гаусса-Зейделя (w = 1)
vector<double> GaussSeidelMethod(const vector<vector<double>>& A, const vector<double>& b, int n) {
    vector<double> X(n, 0.0); // Начальное приближение

    int k = 0;
    while (k < kmax) {
        for (int i = 0; i < n; i++) {
            double sum = SumAx(A, X, i);
            X[i] = (b[i] - sum) / A[i][i];
        }

        double maxDiff = 0.0;
        for (int i = 0; i < n; i++) {
            maxDiff = max(maxDiff, abs(X[i]));
        }

        k++;

        if (maxDiff < eps) {
            cout << "Gauss-Seidel method: the required accuracy per iteration is achieved on iter " << k << endl;
            return X;
        }
    }

    cout << "Gauss-Seidel method: the maximum number of iterations reached (" << kmax << ")" << endl;
    return X;
}

// Метод релаксации
vector<double> RelaxationMethod(const vector<vector<double>>& A, const vector<double>& b, double w, int n) {
    vector<double> X(n, 0.0); // Начальное приближение

    int k = 0;
    while (k < kmax) {
        vector<double> X_new(n);
        for (int i = 0; i < n; i++) {
            double sum1 = SumAx(A, X_new, i);
            double sum2 = SumAx(A, X, i);
            X_new[i] = X[i] + w * ((b[i] - sum1) / A[i][i] - X[i] + sum2);
        }


        double maxDiff = 0.0;
        for (int i = 0; i < n; i++) {
            maxDiff = max(maxDiff, abs(X_new[i] - X[i]));
        }

        X = X_new;
        k++;

        if (maxDiff < eps) {
            cout << "Relaxation method (w = " << w << "): the required accuracy per iteration has been achieved on iter " << k << endl;
            return X;
        }
    }

    cout << "Relaxation method (w = " << w << "): the maximum number of iterations has been reached (" << kmax << ")" << endl;
    return X;
}

void printMatrix(const vector<vector<double>>& matrix, const vector<double>& b, int n) {
    int max_rows = 8;
    int max_cols = 8;

    cout << fixed << setprecision(6);

    for (int i = 0; i < min(n, max_rows); i++) {
        for (int j = 0; j < min(n, max_cols); j++) {
            cout << setw(5) << matrix[i][j] << " ";
        }

        if (n > max_cols) {
            cout << " ...";
        }

        cout << " |     " << b[i] << endl;
    }

    if (n > max_rows) {
        cout << setw(5) << "...";
        for (int j = 1; j < max_cols; j++) {
            cout << " ...";
        }

        cout << " ...";

        cout << " | " << setw(7) << "..." << endl;
    }

    cout << endl;
    cout << fixed << setprecision(8);
}

// Функция для вычисления кубической нормы разности между двумя векторами
double CubicNorm(const vector<double>& x1, const vector<double>& x2) {
    int n = x1.size();
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        norm += pow(abs(x1[i] - x2[i]), 3);
    }
    return cbrt(norm);
}

// Функция для вычисления относительной погрешности
double RelativeError(const vector<double>& x_true, const vector<double>& x_approx) {
    int n = x_true.size();
    double max_abs_x_true = 0.0;
    for (int i = 0; i < n; i++) {
        max_abs_x_true = max(max_abs_x_true, abs(x_true[i]));
    }

    double error = 0.0;
    for (int i = 0; i < n; i++) {
        error = max(error, abs(x_true[i] - x_approx[i]) / max_abs_x_true);
    }
    return error;
}

std::vector<std::vector<double>> GenerateModifiedMatrixA(int n) {
    std::vector<std::vector<double>> Am(n, std::vector<double>(n, 0.0));

    // Строгое диагональное преобладание в первой строке
    for (int j = 0; j < n; j++) {
        Am[0][j] = 0.001 / std::pow(j + 1, 0.3);
    }
    Am[0][0] = 10.0 * std::sqrt(1 + 3);

    // Равенство abs(Aii) = sum(Aij) в остальных строках
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Am[i][j] = 0.001 / std::pow(j + 1, 0.3);
        }
    }

    return Am;
}

int main() {
    // Задайте матрицу A и вектор b
    vector<vector<double>> A = GenerateMatrixA(N);
    vector<double> b = GenerateVectorB(A, N);
    vector<vector<double>> Am = GenerateModifiedMatrixA(N);

    // Вывод матрицы A и вектора b
    cout << "Matrix A:" << endl;
    printMatrix(A, b, N);

    bool isStrictlyDiagonallyDominant = IsStrictlyDiagonallyDominant(A);

    if (isStrictlyDiagonallyDominant) {
        cout << "The matrix is strictly diagonally dominant." << endl;
    }
    else {
        cout << "The matrix is not strictly diagonally dominant." << endl;
    }
    cout << endl << endl;

    // Решение системы методами
    vector<double> solutionJacobi = JacobiMethod(A, b, N);
    vector<double> solutionGaussSeidel = GaussSeidelMethod(A, b, N);
    vector<double> solutionRelaxation1 = RelaxationMethod(A, b, 0.5, N);
    vector<double> solutionRelaxation2 = RelaxationMethod(A, b, 1.5, N);

    // Вычисление точного решения (x*)
    vector<double> x_true(b);

    // Вывод точного решения
    cout << "Precise solution (x*): ";
    cout << fixed << setprecision(6);
    for (double x : x_true) {
        cout << x << " ";
    }
    cout << endl << endl << endl;

    // Вывод приближенных решений
    cout << "Approximate solution by the Jacobi method: ";
    cout << fixed << setprecision(6);
    for (double x : solutionJacobi) {
        cout << x << " ";
    }
    cout << endl;


    cout << "Approximate solution by the Gauss-Seidel method: ";
    cout << fixed << setprecision(6);
    for (double x : solutionGaussSeidel) {
        cout << x << " ";
    }
    cout << endl;

    cout << "Approximate solution by the relaxation method (w = 0.5): ";
    cout << fixed << setprecision(6);
    for (double x : solutionRelaxation1) {
        cout << x << " ";
    }
    cout << endl;

    cout << "Approximate solution by the relaxation method (w = 1.5): ";
    cout << fixed << setprecision(6);
    for (double x : solutionRelaxation2) {
        cout << x << " ";
    }
    cout << endl << endl << endl;

    // Вычисление кубической нормы разности и относительной погрешности
    double normJacobi = CubicNorm(x_true, solutionJacobi);
    double normGaussSeidel = CubicNorm(x_true, solutionGaussSeidel);
    double normRelaxation1 = CubicNorm(x_true, solutionRelaxation1);
    double normRelaxation2 = CubicNorm(x_true, solutionRelaxation2);

    double errorJacobi = RelativeError(x_true, solutionJacobi);
    double errorGaussSeidel = RelativeError(x_true, solutionGaussSeidel);
    double errorRelaxation1 = RelativeError(x_true, solutionRelaxation1);
    double errorRelaxation2 = RelativeError(x_true, solutionRelaxation2);

    // Вывод результатов
    cout << "Cubic norm of difference (Jacobi method): " << normJacobi << endl;
    cout << "Relative error (Jacobi method): " << errorJacobi << endl;

    cout << "Cubic norm of difference (Gauss-Seidel method): " << normGaussSeidel << endl;
    cout << "Relative error (Gauss-Seidel method): " << errorGaussSeidel << endl;

    cout << "Cubic norm of difference (relaxation method, w = 0.5): " << normRelaxation1 << endl;
    cout << "Relative error (relaxation method, w = 0.5): " << errorRelaxation1 << endl;

    cout << "Cubic norm of difference (relaxation method, w = 1.5): " << normRelaxation2 << endl;
    cout << "Relative error (relaxation method, w = 1.5): " << errorRelaxation2 << endl;

    return 0;
}
