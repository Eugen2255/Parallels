#include <iostream>
#define _USE_MATH_DEFINES 
#include <cmath>
#include <vector>


using type =
    #ifdef DOUBLE
        double
    #else
        float
    #endif
    ;

int main() {
    const int size = 10000000;

    std::vector<type> vec (size);

    type sum = 0;

    #ifdef DOUBLE
        double step = 2.0 * M_PI / size;

        for (size_t i = 0; i < size; ++i) {
            vec[i] = sin(step * i);
            sum += vec[i];
        }
        std::cout << "Double: ";

    #else
        float step = 2.0f * M_PI / size;

        for (size_t i = 0; i != size; ++i) {
            vec[i] = sinf(step * i);
            sum += vec[i];
        }
        std::cout << "Float: ";
    #endif
    
    std::cout << "sum = " << sum << std::endl;
    return 0;
}