#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>

bool verify_file(const std::string& filename, int type) {
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "❌ Ошибка: не удалось открыть файл " << filename << "\n";
        return false;
    }

    std::string line;
    size_t lines = 0;
    const double EPS = 1e-3; // Допустимая погрешность для float

    while (std::getline(in, line)) {
        std::istringstream iss(line);
        size_t id;
        double arg1, arg2 = 0.0, res, expected;
        iss >> id >> arg1;
        if (type == 3) iss >> arg2; // pow имеет два аргумента
        iss >> res;

        switch (type) {
            case 1: expected = std::sin(arg1); break;
            case 2: expected = std::sqrt(arg1); break;
            case 3: expected = std::pow(arg1, arg2); break;
        }

        if (std::abs(res - expected) > EPS) {
            std::cerr << "❌ Расхождение в " << filename 
                      << " (ID=" << id << "): ожидалось " << expected 
                      << ", получено " << res << "\n";
            return false;
        }
        ++lines;
    }
    std::cout << "✅ Тест " << filename << " пройден (" << lines << " записей).\n";
    return true;
}

int main() {
    std::cout << "Запуск теста проверки результатов...\n";
    bool ok = true;
    ok &= verify_file("results_sin.txt", 1);
    ok &= verify_file("results_sqrt.txt", 2);
    ok &= verify_file("results_pow.txt", 3);
    
    if (ok) std::cout << "\n🎉 Все тесты успешно пройдены!\n";
    else    std::cout << "\n💥 Тесты завершены с ошибками.\n";
    
    return ok ? 0 : 1;
}