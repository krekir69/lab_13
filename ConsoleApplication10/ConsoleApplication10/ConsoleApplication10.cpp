#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;
const double EPS = 0.0000000000000001;
bool gaussJordan1(vector<vector<double>>& A, vector<double>& X, int m, int n) {
    for (int i = 0; i < n; ++i) {
        // Поиск главной строки
        int v = i;
        for (int j = i + 1; j < m; ++j)
            if (abs(A[j][i]) > abs(A[v][i]))
                v = j;

        if (abs(A[v][i]) < EPS)
            return false; // Нет единственного решения

        // Перестановка строк
        swap(A[i], A[v]);

        // превращаем число в 1
        double div = A[i][i];
        for (int j = 0; j <= n; ++j)
            A[i][j] /= div;

        // Обнуление других строк
        for (int j = 0; j < m; ++j) {
            if (j != i) {
                double c = A[j][i];
                for (int k = 0; k <= n; ++k)
                    A[j][k] -= c * A[i][k];
            }
        }
    }

    for (int i = 0; i < n; ++i)
        X[i] = A[i][n];

    return true;
}
bool gaussJordan2(vector<vector<double>>& A, vector<double>& X, int m, int n) {
    srand(time(0));
    int row = 0;
    for (int i = 0; i < n; ++i) {
        // Поиск главной строки
        int v = row;
        for (int j = i + 1; j < m; ++j)
            if (abs(A[j][i]) > abs(A[v][i]))
                v = j;

        // Перестановка строк
        swap(A[i], A[v]);
        // Обнуление нижних других строк
        for (int j = row + 1; j < m; ++j) {
            double c = A[j][i] / (A[row][i] + EPS);
            for (int k = i; k <= n; ++k)
                A[j][k] -= c * (A[row][k]);

        }
        row++;
    }
    for (int i = n - 1; i >= 0; --i) {
        // Если диагональный элемент ноль — свободная переменная
        if (abs(A[i][i]) < EPS) {
            cout << "\nПеременная x" << i + 1 << " свободная.";
            X[i] = rand() % 10;
            cout << "\nX" << i + 1 << " = " << X[i];

            continue;
        }
        double sum = A[i][n];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[i][j] * X[j];
        }

        X[i] = sum / A[i][i];
    }

    return true;

}
// Печать матрицы
void printMatrix(ofstream& out, const vector<vector<double>>& A, int m, int n) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j <= n; ++j) {
            out << setw(4) << round((A[i][j]) * 100) / 100 << "\t";
        }
        out << "\n";
    }
}
int main() {
    setlocale(LC_ALL, "Ru");
    ifstream inFile("input.txt");
    ofstream outFile("output.txt");

    if (!inFile || !outFile) {
        cout << "Ошибка открытия файлов!\n";
        return 1;
    }

    int m, n,g;
    inFile >> m >> n >> g;

    vector<vector<double>> A(m, vector<double>(n + 1));
    vector<double> X(n);

    if (n > m)
        outFile << "Система несовместна = не имеет решения.\n";
    else
    {
        if (g == 1)
        {
            for (int i = 0; i < m; ++i)
                for (int j = 0; j <= n; ++j)
                    inFile >> A[i][j];

            outFile << "Матрица после приведения к ступенчатому виду:\n";
            bool solved = gaussJordan1(A, X, m, n);
            printMatrix(outFile, A, m, n);
            outFile << "\n";

            if (solved) {
                outFile << "Решение системы:\n";
                for (int i = 0; i < n; ++i)
                    outFile << "x" << (i + 1) << " = " << setw(4) << round((X[i]) * 100.0) / 100.0 << "\n";
            }
            else {
                outFile << "Система не имеет решения.\n";
            }

        }
        else
        {


            for (int i = 0; i < m; ++i)
                for (int j = 0; j <= n; ++j)
                    inFile >> A[i][j];

            outFile << "Матрица после приведения к ступенчатому виду:\n";
            bool solved = gaussJordan2(A, X, m, n);
            printMatrix(outFile, A, m, n);
            outFile << "\n";

            if (solved) {
                outFile << "Решение системы:\n";
                for (int i = 0; i < n; ++i)
                    outFile << "x" << (i + 1) << " = " << setw(4) << round((X[i]) * 100.0) / 100.0 << "\n";
            }
            else {
                outFile << "Система не имеет решения.\n";
            }
        }
    }


    inFile.close();
    outFile.close();

    return 0;
}