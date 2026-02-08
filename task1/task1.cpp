#include <iostream>
#define _USE_MATH_DEFINES 
#include <cmath>
#include <vector>

//#define PI 3,14159265358979323846

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

    type step = 2 * M_PI / size;
    
    type sum = 0;

    #ifdef DOUBLE
    for (size_t i = 0; i < size; ++i) {
        vec[i] = sin(step * i);
        sum += vec[i];
    }
    std::cout << "Double: ";
    #else
    for (size_t i = 0; i != size; ++i) {
        vec[i] = sinf(step * i);
        sum += vec[i];
    }
    std::cout << "Float: ";
    #endif
    

    std::cout << "sum = " << sum << std::endl;
    return 0;
}