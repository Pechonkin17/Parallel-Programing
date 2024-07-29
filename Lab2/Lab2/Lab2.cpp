/*
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;

// Функція для заповнення матриці випадковими числами
void fillMatrix(vector<vector<int>>& matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = rand() % 20 - 10;
        }
    }
}

// Функція для виведення матриці
void printMatrix(const vector<vector<int>>& matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
}

// Функція для підрахунку кількості від'ємних елементів у кожному рядку матриці
vector<int> countNegativeElements(const vector<vector<int>>& matrix, int n) {
    vector<int> negativeCount(n, 0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (matrix[i][j] < 0) {
                negativeCount[i]++;
            }
        }
    }

    return negativeCount;
}

int main() {
    srand(time(0)); // Ініціалізуємо генератор випадкових чисел з часом

    int n = 8; //   <--------------------------------- CHANGE N, and COUNT OF TREADS

    // Створюємо квадратні матриці А та В розміром n*n
    vector<vector<int>> matrixA(n, vector<int>(n));
    vector<vector<int>> matrixB(n, vector<int>(n));

    // Заповнюємо матриці випадковими числами
    fillMatrix(matrixA, n);
    fillMatrix(matrixB, n);

    // Виводимо матриці A
    cout << "Matrix A:" << endl;
    printMatrix(matrixA, n);
    cout << endl;

    // Виводимо матриці B
    cout << "Matrix B:" << endl;
    printMatrix(matrixB, n);
    cout << endl;

    // Підраховуємо кількість від'ємних елементів у кожному рядку матриці A
    vector<int> negativeCountA = countNegativeElements(matrixA, n);
    cout << "Negative numbers in the row, matrix A:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "row " << i + 1 << " have: " << negativeCountA[i] << " negative numbers" << endl;
    }
    cout << endl;

    // Підраховуємо кількість від'ємних елементів у кожному рядку матриці B
    vector<int> negativeCountB = countNegativeElements(matrixB, n);
    cout << "Negative numbers in the row, matrix B:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "row " << i + 1 << " have: " << negativeCountB[i] << " negative numbers" << endl;
    }
    cout << endl;

    return 0;
}
*/















#include <iostream>
#include <vector>
#include <omp.h>
using namespace std;

// Функція для обчислення часу в секундах
double get_time() {
    return omp_get_wtime();
}

// Функція для заповнення матриці випадковими числами
void fillMatrix(vector<vector<int>>& matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = rand() % 20 - 10;
        }
    }
}

// Функція для виведення матриці
void printMatrix(const vector<vector<int>>& matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
}

// Функція для підрахунку кількості від'ємних елементів у кожному рядку матриці A
vector<int> countNegativeElementsMatrixStatically(const vector<vector<int>>& matrix, int n) {
    vector<int> negativeCount(n, 0);
    #pragma omp parallel for schedule(static, 6)
    for (int i = 0; i < n; ++i) {
        int localCount = 0;

        for (int j = 0; j < n; ++j) {
            if (matrix[i][j] < 0) {
                localCount++;
            }
            // cout << "Thread " << omp_get_thread_num() << " (static, 6)\n";
        }

        #pragma omp atomic
        negativeCount[i] += localCount;
    }

    return negativeCount;
}

// Функція для підрахунку кількості від'ємних елементів у кожному рядку матриці B
vector<int> countNegativeElementsMatrixDynamically(const vector<vector<int>>& matrix, int n) {

    vector<int> negativeCount(n, 0);

    #pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < n; ++i) {
        int localCount = 0;

        for (int j = 0; j < n; ++j) {
            if (matrix[i][j] < 0) {
                localCount++;
            }
            // cout << "Thread " << omp_get_thread_num() << " (dynamic, 4)\n";
        }

        #pragma omp atomic
        negativeCount[i] += localCount;
    }

    return negativeCount;

}

// Функція для підрахунку кількості від'ємних елементів у кожному рядку матриці C
vector<int> countNegativeElementsMatrixSequentially(const vector<vector<int>>& matrix, int n) {

    vector<int> negativeCount(n, 0);
    for (int i = 0; i < n; ++i) {

        for (int j = 0; j < n; ++j) {
            if (matrix[i][j] < 0) {
                negativeCount[i]++;
            }
        }
    }

    return negativeCount;

}

int main() {
    srand(time(0)); // Ініціалізуємо генератор випадкових чисел з часом

    int n = 6000; //   <--------------------------------- CHANGE N, and COUNT OF TREADS
    omp_set_num_threads(n);

    // Створюємо квадратні матриці А та В розміром n*n
    vector<vector<int>> matrixA(n, vector<int>(n));
    vector<vector<int>> matrixB(n, vector<int>(n));

    // Заповнюємо матриці випадковими числами
    fillMatrix(matrixA, n);
    fillMatrix(matrixB, n);

  
    // Виводимо матриці A
    cout << "Matrix A:" << endl;
    //printMatrix(matrixA, n);
    cout << endl;

    // Виводимо матриці B
    cout << "Matrix B:" << endl;
    //printMatrix(matrixB, n);
    cout << endl;
    

    // Підраховуємо кількість від'ємних елементів у кожному рядку матриці A і рахуємо час виконання програми
    double start_time = get_time();
    vector<int> C = countNegativeElementsMatrixStatically(matrixA, n);
    double end_time = get_time();
    double static_time = end_time - start_time;
    cout << "Static Time: " << static_time << " seconds" << endl;
    
    cout << "Negative numbers in the row, matrix A, solved with static 6:" << endl;
    for (int i = 0; i < n; ++i) {
       // cout << "row " << i + 1 << " have: " << C[i] << " negative numbers" << endl;
    }
    cout << endl;
 

    // Підраховуємо кількість від'ємних елементів у кожному рядку матриці B і рахуємо час виконання програми
    start_time = get_time();
    vector<int> D = countNegativeElementsMatrixDynamically(matrixB, n);
    end_time = get_time();
    double dynamic_time = end_time - start_time;
    cout << "Dynamic Time: " << dynamic_time << " seconds" << endl;
 
    cout << "Negative numbers in the row, matrix B, solved with dynamic 4:" << endl;
    for (int i = 0; i < n; ++i) {
       // cout << "row " << i + 1 << " have: " << D[i] << " negative numbers" << endl;
    }
    cout << endl;


    // Підраховуємо кількість від'ємних елементів у кожному рядку матриці C і рахуємо час виконання програми
    start_time = get_time();
    vector<int> E = countNegativeElementsMatrixSequentially(matrixA, n);
    end_time = get_time();
    double sequential_time_1 = end_time - start_time;
    cout << "Sequential Time: " << sequential_time_1 << " seconds" << endl;

  
    cout << "Negative numbers in the row, matrix C:" << endl;
    for (int i = 0; i < n; ++i) {
       // cout << "row " << i + 1 << " have: " << E[i] << " negative numbers" << endl;
    }
    cout << endl;
  

    // Підраховуємо кількість від'ємних елементів у кожному рядку матриці D і рахуємо час виконання програми
    start_time = get_time();
    vector<int> F = countNegativeElementsMatrixSequentially(matrixB, n);
    end_time = get_time();
    double sequential_time_2 = end_time - start_time;
    cout << "Sequential Time: " << sequential_time_2 << " seconds" << endl;


    cout << "Negative numbers in the row, matrix D:" << endl;
    for (int i = 0; i < n; ++i) {
      //  cout << "row " << i + 1 << " have: " << F[i] << " negative numbers" << endl;
    }
    cout << endl;
    

    double speedUpA = sequential_time_1 / static_time;
    double speedUpB = sequential_time_2 / dynamic_time;
    double efficiencyA = speedUpA / 4;
    double efficiencyB = speedUpB / 4;

    cout << "speed up A:" << speedUpA << endl;
    cout << "speed up B:" << speedUpB << endl;
    cout << "efficiency A:" << efficiencyA << endl;
    cout << "efficiency B:" << efficiencyB << endl;

    return 0;
}
