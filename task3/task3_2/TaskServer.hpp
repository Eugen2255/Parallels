#pragma once
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <unordered_map>
#include <functional>
#include <atomic>
#include <stdexcept>

template<typename T>
class TaskServer {
public:
    // 1) Запуск сервера
        void start() {
        running.store(true);
        worker = std::jthread([this](std::stop_token stoken) {
            server_loop(stoken);
        });
    }

    // 2) Остановка сервера — ИСПРАВЛЕНО
    void stop() {
        if (running.exchange(false)) {
            // 1. Сначала запрашиваем остановку у jthread
            worker.request_stop();
            // 2. Пробуждаем поток, чтобы он увидел stop_requested()
            queue_cv.notify_all();
            // 3. Ждём завершения
            worker.join();
        }
    }

    // 3) Добавление задачи. Принимает любой callable, возвращающий T
    template<typename F>
    size_t add_task(F&& func) {
        size_t id = next_id.fetch_add(1);
        
        {
            std::lock_guard lock(queue_mutex);
            // Храним как std::function для простоты и универсальности
            task_queue.emplace(id, std::function<T()>(std::forward<F>(func)));
        }
        queue_cv.notify_one();
        return id;
    }

    // 4) Получение результата (блокирующий)
    T request_result(size_t id) {
        std::unique_lock lock(result_mutex);
        result_cv.wait(lock, [this, id] { 
            return results.count(id) > 0 || !running.load(); 
        });
        
        auto it = results.find(id);
        if (it == results.end()) {
            throw std::runtime_error("Result not found for id: " + std::to_string(id));
        }
        
        T res = std::move(it->second);
        results.erase(it); // Удаляем после получения
        return res;
    }

    ~TaskServer() { stop(); }

private:
    void server_loop(std::stop_token stoken) {
        while (!stoken.stop_requested()) {
            std::pair<size_t, std::function<T()>> task_item;
            bool has_task = false;
            
            {
                std::unique_lock lock(queue_mutex);
                queue_cv.wait(lock, [&stoken, this] {
                    return !task_queue.empty() || stoken.stop_requested();
                });
                
                if (stoken.stop_requested() && task_queue.empty()) 
                    break;
                
                if (!task_queue.empty()) {
                    task_item = std::move(task_queue.front());
                    task_queue.pop();
                    has_task = true;
                }
            }
            
            if (!has_task) 
                continue;
            
            // Выполняем задачу БЕЗ захвата мьютекса очереди
            T result = task_item.second();
            
            // Сохраняем результат
            {
                std::lock_guard r_lock(result_mutex);
                results[task_item.first] = std::move(result);
            }
            result_cv.notify_all();
        }
    }

    // Очередь задач: (ID, функция)
    std::queue<std::pair<size_t, std::function<T()>>> task_queue;
    // Контейнер результатов: O(1) вставка/поиск/удаление
    std::unordered_map<size_t, T> results;
    
    std::mutex queue_mutex;
    std::mutex result_mutex;
    std::condition_variable queue_cv;
    std::condition_variable result_cv;
    
    std::jthread worker;
    std::atomic<bool> running{false};
    std::atomic<size_t> next_id{0};
};