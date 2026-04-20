#include "TaskServer.hpp"
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <thread>
#include <string>
#include <chrono>


template<typename T> T fun_sin(T arg) { return std::sin(arg); }
template<typename T> T fun_sqrt(T arg) { return std::sqrt(arg); }
template<typename T> T fun_pow(T x, T y) { return std::pow(x, y); }

void client_sin(TaskServer<double>& srv, size_t N, const std::string& file) {
    std::ofstream out(file);
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dis(-10.0, 10.0);

    for (size_t i = 0; i < N; ++i) {
        double arg = dis(gen);
        // Создаём задачу, замыкающая захватывает аргумент
        auto task = [arg]() { return fun_sin(arg); };
        size_t id = srv.add_task(std::move(task));
        double res = srv.request_result(id);
        out << id << " " << arg << " " << res << "\n";
    }
}

void client_sqrt(TaskServer<double>& srv, size_t N, const std::string& file) {
    std::ofstream out(file);
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dis(0.0, 100.0);

    for (size_t i = 0; i < N; ++i) {
        double arg = dis(gen);
        auto task = [arg]() { return fun_sqrt(arg); };
        size_t id = srv.add_task(std::move(task));
        double res = srv.request_result(id);
        out << id << " " << arg << " " << res << "\n";
    }
}

void client_pow(TaskServer<double>& srv, size_t N, const std::string& file) {
    std::ofstream out(file);
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dis_base(0.1, 5.0);
    std::uniform_real_distribution<double> dis_exp(0.0, 3.0);

    for (size_t i = 0; i < N; ++i) {
        double x = dis_base(gen);
        double y = dis_exp(gen);
        auto task = [x, y]() { return fun_pow(x, y); };
        size_t id = srv.add_task(std::move(task));
        double res = srv.request_result(id);
        out << id << " " << x << " " << y << " " << res << "\n";
    }
}

int main() {
    const size_t N = 100000; 
    std::cout << "Starting server...\n";
    
    auto start = std::chrono::high_resolution_clock::now();
  
    TaskServer<double> server;
    server.start();

    std::thread c1(client_sin, std::ref(server), N, "results_sin.txt");
    std::thread c2(client_sqrt, std::ref(server), N, "results_sqrt.txt");
    std::thread c3(client_pow, std::ref(server), N, "results_pow.txt");

    c1.join();
    c2.join();
    c3.join();

    std::cout << "All clients finished. Stopping server...\n";
    server.stop();
    std::cout << "Server stopped. Files generated successfully.\n";

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Total time: " 
            << std::chrono::duration<double>(end - start).count() 
            << " ms\n";

    return 0;
}