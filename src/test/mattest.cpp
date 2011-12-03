/*
 * =====================================================================================
 *
 *       Filename:  mattest.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/27/2011 08:41:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */


#include <iostream>
#include <typeinfo>
#include <limits>
#include <cstdlib>
#include <ctime>

#include "../classes/matrix"


int main() {

    typedef int T1;
    typedef short T2;
    typedef double T3;

    using namespace boss14420;

#ifndef NDEBUG

//    std::cout << "RAND_MAX = " << RAND_MAX << std::endl;

    matrix::Matrix<T1> m1 = matrix::Matrix<T1>::Random(3,3,-5,5);
    matrix::Matrix<T2> m2 = matrix::Matrix<T2>::Random(3,3,-5,5);
    matrix::Matrix<T3> m3 = matrix::Matrix<T3>::Random(3,3,-5,5);

    typedef decltype(T1()+(T2()*T3())+T2()) TRes;
//    Matrix<TRes> mt = mi + md + ms;
    
    std::cout << "Matrix<" << typeid(T1).name() << "> m1 = \n" << m1 << std::endl;
    std::cout << "Matrix<" << typeid(T2).name() << "> m2 = \n" << m2 << std::endl;
    std::cout << "Matrix<" << typeid(T3).name() << "> m3 = \n" << m3 << std::endl;
    
    std::cout << "Matrix<" << typeid(TRes).name()
        << "> m1 + (m2 * m3) + m3 = \n" << m1+(m2*m3)+m3 << std::endl;
#endif

    return 0;
}
