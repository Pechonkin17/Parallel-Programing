#include <iostream>
#include <cmath>
#include <omp.h>

using namespace std;

const int n = 2000;

// Функція для обчислення часу в секундах
double get_time() {
    return omp_get_wtime();
}

// Функція для множення матриці на вектор
void matrix_vector_multiply_parallel(double** matrix, double* vector, double* result) {
#pragma omp parallel for 
    for (int i = 0; i < n; ++i) {
        result[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

void matrix_vector_multiply_sequential(double** matrix, double* vector, double* result) {
    for (int i = 0; i < n; ++i) {
        result[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

int main() {
    cout << "n = " << n << endl;

    // Ініціалізація матриці та вектора
    double** matrix = new double* [n];
    double* vector = new double[n];
    double* result = new double[n];

    for (int i = 0; i < n; ++i) {
        matrix[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = (pow(j, 2) + 3 * j + 2) / (pow(i + 1, 1/2));
        }
        vector[i] = log10(pow(i, 1/2) * tan(2 * i)) / (sin(i) + 3);
    }

    // Послідовний спосіб обчислення
    double start_time = get_time();
    matrix_vector_multiply_sequential(matrix, vector, result);
    double end_time = get_time();
    double sequential_time = end_time - start_time;
    cout << "Sequential Time: " << sequential_time << " seconds" << endl;

    // Паралельний спосіб обчислення
    int num_threads = 4; // <---------------------------------------------------------- THREADS CHANGE
    omp_set_num_threads(num_threads);

    start_time = get_time();
    matrix_vector_multiply_parallel(matrix, vector, result);
    end_time = get_time();
    double parallel_time = end_time - start_time;
    cout << "Parallel Time (" << num_threads << " threads): " << parallel_time << " seconds" << endl;

    cout << "Efficiency: " << sequential_time / parallel_time << endl;

    cout << "Speedup: " << (sequential_time / parallel_time) / 4;

    // Очищеня пам'яті
    for (int i = 0; i < n; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] vector;
    delete[] result;

    return 0;
}
