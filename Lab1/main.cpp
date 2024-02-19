#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;

vector <double> tridiagonal_solve(const vector<vector<double>>& A, const vector<double>& d) {
    int n = d.size();
    vector<double> a(n);
    vector<double> b(n);
    vector<double> c(n);

    for (int i = 0; i < n; i++) {
        if (i != 0) a[i] = A[i][i - 1];
        b[i] = A[i][i];
        if (i != n - 1) c[i] = A[i][i + 1];
    }

    vector<double> p(n, 0);
    vector<double> q(n, 0);
    vector<double> x(n, 0);

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        p[i] = -c[i] / (a[i] * p[i - 1] + b[i]);
        q[i] = (d[i] - a[i] * q[i - 1]) / (a[i] * p[i - 1] + b[i]);
    }

    x[n - 1] = q[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = p[i] * x[i + 1] + q[i];
    }
    return x;
}


pair<vector<vector<double>>, vector<vector<double>>> lu_decomposition(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<vector<double>> U(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            U[i][k] = A[i][k];
            for (int j = 0; j < i; j++) {
                U[i][k] -= L[i][j] * U[j][k];
            }
        }
        L[i][i] = 1.0;
        for (int k = i + 1; k < n; k++) {
            L[k][i] = A[k][i] / U[i][i];
            for (int j = 0; j < i; j++) {
                L[k][i] -= L[k][j] * U[j][i] / U[i][i];
            }
        }
    }
    return { L, U };
}
vector<vector<double>>  inverse_matrix(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> inv(n, vector<double>(n, 0.0));
    pair<vector<vector<double>>, vector<vector<double>>> L, U = lu_decomposition(A);
    for (int i = 0; i < n; i++) {
        vector<double> e(n, 0.0);
        e[i] = 1.0;
        vector<double> y(n, 0.0);
        for (int j = 0; j < n; j++) {
            y[j] = e[j];
            for (int k = 0; k < j; k++) {
                y[j] -= lu_decomposition(A).first[j][k] * y[k];
            }
        }
        vector<double> x = tridiagonal_solve(lu_decomposition(A).second, y);  // Используется метод прогонки
        for (int j = 0; j < n; j++) {
            inv[j][i] = x[j];
        }
    }
    return inv;
}
double condition_number(const vector<vector<double>>& A, const vector<vector<double>>& A_inv) {
    double norm_A = 0;
    double norm_A_inv = 0;

    for (int i = 0; i < A.size(); i++) {
        double row_sum_A = 0.0;
        double row_sum_A_inv = 0.0;
        for (int j = 0; j < A.size(); j++) {
            row_sum_A += abs(A[i][j]);
            row_sum_A_inv += abs(A_inv[i][j]);
        }
        norm_A = max(norm_A, row_sum_A);
        norm_A_inv = max(norm_A_inv, row_sum_A_inv);
    }

    return norm_A * norm_A_inv;
}



int main() {
    int n = 15;
    cout << "Matrix:";
    cout << endl;
    vector<vector<double>> A(n, vector<double>(n, 0));
    // Заполнение главной диагонали
    for (int i = 0; i < n; i++) {
        A[i][i] = 5 * (i + 1);
    }
    // Заполнение диагонали над главной
    for (int i = 0; i < n-1; i++) {
        A[i][i + 1] = -((i + 1) * sqrt((i + 1) + 1));
    }
    // Заполнение диагонали под главной
    for (int i = 1; i < n; i++) {
        A[i][i - 1] = -((i + 1) * sqrt((i + 1) - 1));
    }
    // Вывод матрицы
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setprecision(5) << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Matrix d:";
    cout << endl;
    //Заполнение и вывод матрицы В
    vector <double> d(n);
    for (int i = 0; i < n; i++)
    {
        d[i] = 3 * sqrt((i + 1));
        cout << setprecision(5) << d[i] << endl;
    }
    cout << endl;
    cout << "Run-through method:" << endl;
    vector <double> result = tridiagonal_solve(A, d);
    for (int i = 0; i < n-1; i++) {
        cout << "x" << i + 1 << " = " << result[i] << endl;
    }
    tridiagonal_solve(A, d);
    cout << endl << "LU-decomposition:" << endl;
    pair<vector<vector<double>>, vector<vector<double>>> result_LU = lu_decomposition(A);
    cout << "L:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << result_LU.first[i][j] << " ";
        }
        cout << endl;
    }
    cout << "U:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << result_LU.second[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << "Inverse matrix:" << endl;
    vector <vector <double>> result_2 = inverse_matrix(A);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << result_2[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    cout << "The number of conditionality:" << endl;
    double result_3 = condition_number(A, result_2);
    cout << result_3 << endl << endl;




    ///
  vector<vector<double>> Res(n, vector<double>(n + 1, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++)
            {
                Res[i][j] += result_LU.first[i][k] * result_LU.second[k][j];
            }

        }
    }
    cout << "LU: " << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << setw(4) << Res[i][j];
        }
        cout << endl;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Res[i][j] -= A[i][j];
        }
    }
    double norm = 0.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            norm += Res[i][j] * Res[i][j];
        }
    }
    
    norm = sqrt(norm);
    cout << "Norm of testing: " << norm << endl;
    ///
    



    cout << endl;

    auto start = chrono::high_resolution_clock::now();
    vector<double> x = tridiagonal_solve(A, d);
    auto end = chrono::high_resolution_clock::now();
    double duration = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1e9;
    cout << "Execution time of tridiagonal_solve: " << duration << " seconds\n";

    auto start_1 = chrono::high_resolution_clock::now();
    pair<vector<vector<double>>, vector<vector<double>>> x_1 = lu_decomposition(A);
    auto end_1 = chrono::high_resolution_clock::now();
    double duration_1 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1e9;
    cout << "Execution time of lu_decompositio: " << duration_1 << " seconds\n";

    auto start_2 = chrono::high_resolution_clock::now();
    vector<vector<double>> x_2 = inverse_matrix(A);
    auto end_2 = chrono::high_resolution_clock::now();
    double duration_2 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1e9;
    cout << "Execution time of inverse_matrix: " << duration_2 << " seconds\n";

    auto start_3 = chrono::high_resolution_clock::now();
    double x_3 = condition_number(A, result_2);
    auto end_3 = chrono::high_resolution_clock::now();
    double duration_3 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1e9;
    cout << "Execution time of condition_number: " << duration_3 << " seconds\n";


    return 0;
}
