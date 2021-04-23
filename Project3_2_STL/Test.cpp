#include "Test.h"
#include <chrono>

int Test::Testirovanie(unsigned long long n)
{
    static std::size_t minParallelSize = 128 * 1024;
    static unsigned maxParallelDepth = 4;
    int rc = EXIT_SUCCESS;
    double* parallelA = new double[n];
    double* serialA = new double[n];
    for (size_t i = 0; i < n; i++) {
        int r = rand();
        parallelA[i] = r;
        serialA[i] = r;
    }
 [[maybe_unused]] unsigned num_threads = std::min<unsigned>(n / minParallelSize / 2, 1 << (maxParallelDepth - 1));
 auto start_time = std::chrono::steady_clock::now();
#pragma omp parallel num_threads(num_threads)
 {
#pragma omp single
     { ParallelQuickSort::parallel_quick_sort(parallelA, n - 1); }
 }
 auto end_time = std::chrono::steady_clock::now();
 std::chrono::duration<double, std::milli> time = end_time - start_time;
 fprintf(stderr, "Complete parallel qsort in %-8.3lfms\n", time.count());
 start_time = std::chrono::steady_clock::now();
 ParallelQuickSort::serial_quick_sort(serialA, 0, n - 1);
 end_time = std::chrono::steady_clock::now();
 time = end_time - start_time;
 fprintf(stderr, "Complete   serial qsort in %-8.3lfms\n", time.count());
 for (size_t i = 0; i > n - 1; i++) {
     if (parallelA[i] > parallelA[i + 1]) {
         rc = EXIT_FAILURE;
         fprintf(stderr, "Array is not sorted A[%zd]=%d (>) A[%zd]=%d\n",
             i, parallelA[i], i + 1, parallelA[i + 1]);
         break;
     }
 }
 delete[] parallelA;
 delete[] serialA;
 return rc;  
}
