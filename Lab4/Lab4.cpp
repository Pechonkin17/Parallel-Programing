#include <iostream>
#include <cmath>
#include <mpi.h>
#include <windows.h>

using namespace std;

double f(double x) {
    double a = 18;
    double b = 18 * 2;
    double k = 18 / 2.0;
    double z = 18 * 18;
    return exp(x) - ((pow(z, 2) + 12 * x * z - 3 * pow(x, 2)) / (18 * z - 1));
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = 18 * 10000000;
    double h = 18.0 / N;
    double local_sum = 0.0;
    double global_sum = 0.0;

    DWORD dwStart = GetTickCount();

    for (int i = rank + 1; i <= N + 1; i += size) {
        double x = i * h;
        local_sum += f(x);
    }

    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        DWORD dwRunTime = GetTickCount() - dwStart;
        cout << "Results:" << endl;
        cout << "Total sum: " << global_sum << endl;
        cout << "Time taken: " << dwRunTime << " milliseconds" << endl;
    }

    MPI_Finalize();
    return 0;
}


// mpiexec -n 2 Lab4.exe
// mpiexec -n 4 Lab4.exe