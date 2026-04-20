#include <iostream>
#include <vector>
#include <chrono>
#include <future>

#ifndef M
    #define M 20000
#endif
int m = M;
int n = m;

void compute_chunk(double *a, double *b, double *c, int m, int n, int thread_id, int num_threads) {
    for (int i = thread_id; i < m; i += num_threads) {
        c[i] = 0.0;
        for (int j = 0; j < n; ++j)
            c[i] += a[i * n + j] * b[j];
    }
}

void init_chunk(double *a, double *c, int m, int n, int thread_id, int num_threads) {
    for (int i = thread_id; i < m; i += num_threads) {
        for (int j = 0; j < n; ++j)
            a[i * n + j] = i + j;
        c[i] = 0.0;
    }
}

double run_parallel() {
    std::vector<double> a(m * n), b(n), c(m);
    
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (const char* env_threads = std::getenv("NTHREADS")) {
        int parsed = std::atoi(env_threads);
        if (parsed > 0) num_threads = static_cast<unsigned int>(parsed);
    }
    if (num_threads == 0) num_threads = 4;
        
    std::vector<std::future<void>> init_futures;
    for (unsigned int t = 0; t < num_threads; ++t)
        init_futures.push_back(std::async(std::launch::async, 
                                          init_chunk, a.data(), c.data(), m, n, t, num_threads));
    for (auto &f : init_futures) f.wait();
    
    for (int j = 0; j < n; ++j) b[j] = j;
    
    auto t_start = std::chrono::high_resolution_clock::now();
    std::vector<std::future<void>> compute_futures;
    for (unsigned int t = 0; t < num_threads; ++t)
        compute_futures.push_back(std::async(std::launch::async,
                                             compute_chunk, a.data(), b.data(), c.data(), m, n, t, num_threads));
    for (auto &f : compute_futures) f.wait();
    
    auto t_end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t_end - t_start).count();
}

int main() {
    int ntest = 40;
    double tparallel = 0.0;
    for (int i = 0; i < ntest; ++i)
        tparallel += run_parallel();
    
    std::cout << "Avg time (std::async): " << tparallel / ntest << " sec." << std::endl;
    return 0;
}