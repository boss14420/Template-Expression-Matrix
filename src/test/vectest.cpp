/*
 * =====================================================================================
 *
 *       Filename:  vectest.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/27/2011 04:01:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */

#include "../classes/vector.hpp"
#include <cstdlib>
#include <iostream>

int main() {
    typedef double T;

    T k = 2.0;
    std::vector<T> v1(4), v2(4);

    std::srand(10);
    for(int i = 0; i < 4; ++i) {
        v1[i] = std::rand() % 10 - 5;
        v2[i] = std::rand() % 10 - 5;
    }


    Vec<T> V1(v1), V2(v2);

    auto func = [&]() -> decltype(V1+2.0*V2-3.0*V1) {
        return V1 + 2.0*V2 - 3.0* V1;
    };

    Vec<T> V3 = (V1*V2) * V1 + k * (V2 - V1);

    std::cout << "V1 = (" << V1[0] << ", " << V1[1] << ", " << V1[2] << ", " << V1[3] << ")\n";
    std::cout << "V2 = (" << V2[0] << ", " << V2[1] << ", " << V2[2] << ", " << V2[3] << ")\n";
    std::cout << "V3 = (V1 * V2) * V1 + " << k << " * (V2 - V1) = (" 
        << V3[0] << ", " << V3[1] << ", " << V3[2] << ", " << V3[3] << ")\n";

    auto V4 = func();
    std::cout << "V4 = func() = (" << V4[0] << ", " << V4[1] << ", " << V4[2] << ", " << V4[3] << ")\n";

    auto V5 = V1+V2-V3+V2;
    std::cout << "V5 = V1+V2-V3+V2 = (" << V5[0] << ", " << V5[1] << ", " << V5[2] << ", " << V5[3] << ")\n";

    return 0;
}
