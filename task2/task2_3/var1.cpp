#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <omp.h>

using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

double TAU = 0.000005;
double EPSILON = 1e-5;
int MAX_ITER = 2000;
int N = 20000; 


Mat create_matrix(int n) {
    Mat A(n, Vec(n, 1.0));
    for (int i = 0; i < n; ++i)
        A[i][i] = 2.0;  
    return A;
}

Vec create_vector(int n) {
    return Vec(n, n + 1.0);
}

Vec solve(const Mat& A, const Vec& b, double tau, double eps, int& out_iter) {
    int n = b.size();
    Vec x(n, 0.0);
    
    for (int k = 0; k < MAX_ITER; ++k) {
        Vec x_new(n);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            double Ax_i = 0.0;
            for (int j = 0; j < n; ++j)
                Ax_i += A[i][j] * x[j];
            x_new[i] = x[i] - tau * (Ax_i - b[i]);
        }
        
        double diff_norm = 0.0, x_new_norm = 0.0;
        #pragma omp parallel for reduction(+:diff_norm,x_new_norm)
        for (int i = 0; i < n; ++i) {
            double d = x_new[i] - x[i];
            diff_norm += d * d;
            x_new_norm += x_new[i] * x_new[i];
        }
        
        diff_norm = std::sqrt(diff_norm);
        x_new_norm = std::sqrt(x_new_norm);
        
        if (diff_norm < eps * x_new_norm) {
            out_iter = k + 1;
            return x_new;
        }
        x = std::move(x_new);
    }
    out_iter = MAX_ITER;
    return x;
}

int main() {
    Mat A = create_matrix(N);
    Vec b = create_vector(N);
    
    int iterations = 0;
    int ntest = 30;
    double tparallel = 0;
    for(int i = 0; i < ntest; ++i){
        auto start = std::chrono::high_resolution_clock::now();
        Vec x = solve(A, b, TAU, EPSILON, iterations);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "x0:" << x[0] << "\n";
        //std::cout << "iterations:" << iterations << "\n";
        tparallel += std::chrono::duration<double>(end - start).count();
    }
    
    std::cout << "Avg time: " << tparallel / ntest << std::endl;
    
    
    //std::cout << std::fixed << std::setprecision(6);
    //std::cout << "variant:1\n";
    //std::cout << "time:" << elapsed << "\n";
    //std::cout << "iterations:" << iterations << "\n";
    //std::cout << "x0:" << x[0] << "\n";  
    
    return 0;
}