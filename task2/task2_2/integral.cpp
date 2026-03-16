#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <omp.h>


using Func = std::function<double(double)>;

double func(double x)
{
    return std::exp(-x * x);
}


double integrate_omp(const Func& f, double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.0;

    #pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = n / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? n : (lb + items_per_thread);
        
        double sumloc = 0.0;

        for (int i = lb; i < ub; i++)
            sumloc += f(a + h * (i + 0.5));

        #pragma omp atomic
        sum += sumloc;
    }

    sum *= h;
    return sum;
}

constexpr double PI = 3.14159265358979323846;
constexpr double A = -4.0;
constexpr double B = 4.0;
constexpr int NSTEPS = 80000000;

double get_elapsed_time(std::chrono::high_resolution_clock::time_point start)
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    return diff.count();
}


double run_parallel()
{
    auto start = std::chrono::high_resolution_clock::now();
    double res = integrate_omp(func, A, B, NSTEPS);
    double t = get_elapsed_time(start);
    return t;
}

int main()
{    
    int ntest = 50;
    double tparallel = 0;
    for(int i = 0; i < ntest; ++i)
        tparallel += run_parallel();
    
    std::cout << "Avg time: " << tparallel / ntest << std::endl;
    
    return 0;
}