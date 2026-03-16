#include <iostream>
#include <cstdio>
#include <vector>
#include <chrono>
#include <omp.h>

#ifndef M
    #define M 20000
#endif
int m = M;
int n = m;

/*
* matrix_vector_product: Compute matrix-vector product c[m] = a[m][n] * b[n]
*/
void matrix_vector_product(double *a, double *b, double *c, int m, int n){
    for (int i = 0; i < m; ++i) {
        c[i] = 0.0;
        for (int j = 0; j < n; ++j)
            c[i] += a[i * n + j] * b[j];
    }
}

void run_serial() {
    std::vector<double> a(m * n);
    std::vector<double> b(n);
    std::vector<double> c(m);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            a[i * n + j] = i + j;
    }

    for (int j = 0; j < n; ++j)
        b[j] = j;

    auto t_start = std::chrono::high_resolution_clock::now();
    matrix_vector_product(a.data(), b.data(), c.data(), m, n);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;

    std::cout << "Elapsed time (serial): " << elapsed.count() << " sec." << std::endl;
}

/* matrix_vector_product_omp: Compute matrix-vector product c[m] = a[m][n] * b[n] */
void matrix_vector_product_omp(double *a, double *b, double *c, int m, int n){
    #pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);

        for (int i = lb; i <= ub; ++i) {
            c[i] = 0.0;
            for (int j = 0; j < n; ++j)
                c[i] += a[i * n + j] * b[j];
        }
    }
}

double run_parallel() {
    std::vector<double> a(m * n);
    std::vector<double> b(n);
    std::vector<double> c(m);

    // Параллельная инициализация массивов
    #pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        
        for (int i = lb; i <= ub; ++i) {
            for (int j = 0; j < n; ++j)
                a[i * n + j] = i + j;
            c[i] = 0.0;
        }
    }

    for (int j = 0; j < n; j++)
        b[j] = j;

    auto t_start = std::chrono::high_resolution_clock::now();
    matrix_vector_product_omp(a.data(), b.data(), c.data(), m, n);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;

    //std::cout << "Elapsed time (parallel): " << elapsed.count() << " sec." << std::endl;

    return elapsed.count();
}

int main(int argc, char **argv){
    //std::cout << "Matrix-vector product (c[m] = a[m, n] * b[n]; m = " << m << ", n = " << n << ")\n";

    int ntest = 50;
    double tparallel = 0.0;
    for(int i = 0; i < ntest; ++i){
        //size_t memory_bytes = (static_cast<size_t>(m) * n + m + n) * sizeof(double);
        //size_t memory_mib = memory_bytes / (1024 * 1024);
        //std::cout << "Memory used: " <<   memory_mib << " MiB\n";
        
        //run_serial();
        tparallel += run_parallel();
    }
    std::cout << "Avg time: " << tparallel / ntest << std::endl;

    return 0;
}