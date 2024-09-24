#include <iostream>
#include <pthread.h>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <vector>
#include <thread>  // For getting the number of hardware threads

// Function to calculate the number of digits in an unsigned long long number
int numberOfDigits(unsigned long long n) {
    return n == 0 ? 1 : static_cast<int>(log10(n)) + 1;
}

// Function to compute the rounded square root using Heron's method (adapted for integers)
int roundedSquareRoot(unsigned long long n) {
    int digits = numberOfDigits(n);
    unsigned long long x_k;

    // Initial guess based on the number of digits
    if (digits % 2 == 1) {
        x_k = 2 * pow(10, (digits - 1) / 2);
    } else {
        x_k = 7 * pow(10, (digits - 2) / 2);
    }

    unsigned long long x_k1 = 0;
    int iterations = 0;

    while (true) {
        x_k1 = (x_k + (n + x_k - 1) / x_k) / 2;  // The integer division equivalent of the formula

        iterations++;
        if (x_k == x_k1) break;  // If x_k == x_k1, stop the iteration
        x_k = x_k1;
    }

    return iterations;
}

// Thread arguments struct
struct ThreadArgs {
    unsigned long long startRange;
    unsigned long long endRange;
    unsigned long long totalIterations;
};

// Function that each thread will execute
void* threadFunction(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    unsigned long long totalIterations = 0;

    // Loop over the assigned range and calculate the iterations
    for (unsigned long long n = args->startRange; n < args->endRange; ++n) {
        int iterations = roundedSquareRoot(n);
        totalIterations += iterations;
    }

    args->totalIterations = totalIterations;
    return nullptr;
}

int main() {
    const unsigned long long startRange = 10000000000000ULL;  // Start of the range 10^13
    const unsigned long long endRange = 100000000000000ULL;   // End of the range 10^14 (exclusive)
    unsigned long long numberOfNumbers = endRange - startRange;  // Total numbers in the range

    // Get the number of hardware threads (cores)
    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) {
        numThreads = 4;  // Fallback to 4 threads if hardware_concurrency fails
    }

    // Create an array to hold thread arguments
    std::vector<ThreadArgs> threadArgs(numThreads);
    std::vector<pthread_t> threads(numThreads);

    unsigned long long numbersPerThread = numberOfNumbers / numThreads;

    // Start measuring time
    auto start = std::chrono::high_resolution_clock::now();

    // Create and launch threads
    for (unsigned int i = 0; i < numThreads; ++i) {
        threadArgs[i].startRange = startRange + i * numbersPerThread;
        threadArgs[i].endRange = (i == numThreads - 1) ? endRange : threadArgs[i].startRange + numbersPerThread;
        threadArgs[i].totalIterations = 0;
        pthread_create(&threads[i], nullptr, threadFunction, &threadArgs[i]);
    }

    // Wait for all threads to finish
    unsigned long long totalIterations = 0;
    for (unsigned int i = 0; i < numThreads; ++i) {
        pthread_join(threads[i], nullptr);
        totalIterations += threadArgs[i].totalIterations;
    }

    // Stop measuring time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // Calculate the average number of iterations
    double averageIterations = static_cast<double>(totalIterations) / numberOfNumbers;

    // Output the results
    std::cout << "Average number of iterations for the range [10^13, 10^14): " 
              << std::fixed << std::setprecision(10) << averageIterations << std::endl;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}