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

#include "vector.hpp"
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

    Vec<T> V3 = (V1*V2) * V1 + k * (V2 - V1);

    std::cout << "V1 = (" << V1[0] << ", " << V1[1] << ", " << V1[2] << ", " << V1[3] << ")\n";
    std::cout << "V2 = (" << V2[0] << ", " << V2[1] << ", " << V2[2] << ", " << V2[3] << ")\n";
    std::cout << "V3 = (V1 * V2) * V1 + " << k << " * (V2 - V1) = (" 
        << V3[0] << ", " << V3[1] << ", " << V3[2] << ", " << V3[3] << ")\n";

    return 0;
}
