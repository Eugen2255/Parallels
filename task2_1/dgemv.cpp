#include <iostream>
#include <cstdio>
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
    double* a = new double[m * n];
    double* b = new double[n];
    double* c = new double[m];

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            a[i * n + j] = i + j;
    }

    for (int j = 0; j < n; ++j)
        b[j] = j;

    double t = omp_get_wtime();
    matrix_vector_product(a, b, c, m, n);
    t = omp_get_wtime() - t;

    std::cout << "Elapsed time (serial): " << t << " sec." << std::endl;

    delete[] a;
    delete[] b;
    delete[] c;
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

void run_parallel() {
    double* a = new double[m * n];
    double* b = new double[n];
    double* c = new double[m];

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

    double t = omp_get_wtime();
    matrix_vector_product_omp(a, b, c, m, n);
    t = omp_get_wtime() - t;

    std::cout << "Elapsed time (parallel): " << t << " sec." << std::endl;

    delete[] a;
    delete[] b;
    delete[] c;
}

int main(int argc, char **argv){
    std::cout << "Matrix-vector product (c[m] = a[m, n] * b[n]; m = " << m 
              << ", n = " << n << ")\n";
    
    size_t memory_bytes = (static_cast<size_t>(m) * n + m + n) * sizeof(double);
    size_t memory_mib = memory_bytes / (1024 * 1024);
    
    std::cout << "Memory used: " << memory_mib << " MiB\n";
    
    run_serial();
    run_parallel();

    return 0;
}