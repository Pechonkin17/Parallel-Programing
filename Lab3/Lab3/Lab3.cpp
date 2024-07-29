/*
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;

// Function defining the integrand
double integrand(double x) {
    return 1 / (0.35 * pow(cos(0.6 * x), 2) + pow(sin(0.6 * x), 2));
}

// Function to calculate the integral using left rectangle method
double calculateIntegral(double lowerLimit, double upperLimit, double stepSize) {
    double integral = 0.0;

    // Loop through the intervals with step size and calculate the area of each rectangle
    for (double x = lowerLimit; x < upperLimit; x += stepSize) {
        integral += integrand(x) * stepSize;
    }

    return integral;
}

// Function defining the integrand for the second integral
double integrand2(double x) {
    return (x * sqrt(0.7 * pow(x, 2) + 1.2 * x + 2.1)) / pow((21. + x), 2);
}

// Function to calculate the integral using left rectangle method for the second integral
double calculateIntegral2(double lowerLimit, double upperLimit, double stepSize) {
    double integral = 0.0;

    // Loop through the intervals with step size and calculate the area of each rectangle
    for (double x = lowerLimit; x < upperLimit; x += stepSize) {
        integral += integrand2(x) * stepSize;
    }

    return integral;
}

// Function to write results to a text file
void writeResultsToFile(const string& filename, double lowerLimit, double upperLimit, const double* stepSizes, int numSteps) {
    ofstream outputFile(filename);

    // Writing header to the file
    outputFile << "Lab Work #3: Download and synchronization in OpenMP" << endl;
    outputFile << "Student Group: KN-314" << endl;
    outputFile << "Student's Full Name: Pechonkin Maksym Dmytrovych" << endl;
    outputFile << "Variant Number: 19" << endl;
    outputFile << "Task: Compute integrals" << endl;

    // Performing numerical experiment with different step sizes
    for (int i = 0; i < numSteps; ++i) {
        double stepSize = stepSizes[i];

        outputFile << "The beginning of the closed section with step size " << stepSize << endl;

        double result = calculateIntegral(lowerLimit, upperLimit, stepSize);
        outputFile << "Result of the first integral: " << result << endl;

        double result2 = calculateIntegral2(0.3, 1.3, stepSize);
        outputFile << "Result of the second integral: " << result2 << endl;


        outputFile << "The end of the closed section with step size " << stepSize << endl;
    }

    // Closing the file
    outputFile.close();
}

int main() {
    // Input data
    double lowerLimit = 0.1;
    double upperLimit = 1.1;
    double stepSizes[] = { 0.00002, 0.00005, 0.0001 }; // Different step sizes for discretization
    int numSteps = sizeof(stepSizes) / sizeof(stepSizes[0]);

    int n = 2; //   <--------------------------------- CHANGE N, and COUNT OF TREADS
    omp_set_num_threads(n);


    // Calling function to write results to file
    writeResultsToFile("results.txt", lowerLimit, upperLimit, stepSizes, numSteps);

    return 0;
}
*/









#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;

omp_lock_t lock;

// Function get time
double get_time() {
    return omp_get_wtime();
}

// Function defining the integrand
double integrand(double x) {
    return 1 / (0.35 * pow(cos(0.6 * x), 2) + pow(sin(0.6 * x), 2));
}

// Function to calculate the integral using left rectangle method
double calculateIntegral(double lowerLimit, double upperLimit, double stepSize) {
    double integral = 0.0;
    int numSteps = static_cast<int>((upperLimit - lowerLimit) / stepSize); // Calculate number of steps

    // Parallel loop using OpenMP
    for (int i = 0; i < numSteps; ++i) {
        double x = lowerLimit + i * stepSize; // Convert integer index to floating-point value
        integral += integrand(x) * stepSize;
    }

    return integral;
}

// Function defining the integrand for the second integral
double integrand2(double x) {
    return (x * sqrt(0.7 * pow(x, 2) + 1.2 * x + 2.1)) / pow((21. + x), 2);
}

// Function to calculate the integral using left rectangle method for the second integral
double calculateIntegral2(double lowerLimit, double upperLimit, double stepSize) {
    double integral = 0.0;
    int numSteps = static_cast<int>((upperLimit - lowerLimit) / stepSize); // Calculate number of steps

    // Parallel loop using OpenMP
    for (int i = 0; i < numSteps; ++i) {
        double x = lowerLimit + i * stepSize; // Convert integer index to floating-point value
        integral += integrand2(x) * stepSize;
    }

    return integral;
}

// Function to write results to a text file
void writeResultsToFile(const string& filename, const double* stepSizes, int numSteps) {
    ofstream outputFile(filename);

    // Writing header to the file
    outputFile << "Lab Work #3: Download and synchronization in OpenMP" << endl;
    outputFile << "Student Group: KN-314" << endl;
    outputFile << "Student's Full Name: Pechonkin Maksym Dmytrovych" << endl;
    outputFile << "Variant Number: 18" << endl;
    outputFile << "Task: Compute the integrals" << endl;

    
    int i;
    omp_init_lock(&lock);

    // Performing numerical experiment with different step sizes
    #pragma omp master
    #pragma omp parallel for ordered private (i)
    for (i = 0; i < numSteps; ++i) {
        omp_set_lock(&lock);
        
        //omp_test_lock(&lock);

        double stepSize = stepSizes[i];

        outputFile << "The beginning of the closed section with step size " << stepSize << endl;

        double start_time = get_time();
        double result = calculateIntegral(0.1, 1.1, stepSize);
        double end_time = get_time();
        double static_time = end_time - start_time;
        cout << "time result 1: " << static_time << endl;
        outputFile << "Result of the first integral: " << result << endl;

        start_time = get_time();
        double result2 = calculateIntegral2(0.5, 1.5, stepSize);
        end_time = get_time();
        static_time = end_time - start_time;
        cout << "time result 2: " << static_time << endl;
        outputFile << "Result of the second integral: " << result2 << endl;


        outputFile << "The end of the closed section with step size " << stepSize << endl;

        omp_unset_lock(&lock);
    }
     
    omp_destroy_lock(&lock);

    outputFile.close();
}

int main() {

    // Input data
    double stepSizes[] = { 0.00002, 0.00005, 0.0001 }; // Different step sizes for discretization
    int numSteps = sizeof(stepSizes) / sizeof(stepSizes[0]);

    int n = 2; //   <--------------------------------- CHANGE N, and COUNT OF TREADS
    omp_set_num_threads(n);


    // Calling function to write results to file
    writeResultsToFile("results.txt", stepSizes, numSteps);

    return 0;
}
