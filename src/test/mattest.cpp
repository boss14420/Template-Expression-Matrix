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

#define CREATE_TEST(EXPR) do { std::cout << "----" << #EXPR << "------------------------------------" << std::endl;\
                                std::cout << #EXPR << " = \n" << (EXPR) << std::endl << "-------------" << std::endl;\
                            } while(0);

int main() {

    typedef int T1;
    typedef short T2;
    typedef double T3;
    typedef float T4;

    using namespace boss14420;

//    std::cout << "RAND_MAX = " << RAND_MAX << std::endl;

    matrix::Matrix<T1> m1 = matrix::Matrix<T1>::Random(3,3,-5,5);
    matrix::Matrix<T2> m2 = matrix::Matrix<T2>::Random(3,3,-5,5);
//    matrix::Matrix<T3,matrix::ColumnMajor> m3 = matrix::Matrix<T3,matrix::ColumnMajor>::Random(3,3,-5,5);
    matrix::Matrix<T3> m3 = matrix::Matrix<T3>::Random(3,3,-5,5);
    matrix::Matrix<T4> m4 = matrix::Matrix<T4>::Random(3,3,-5,5);

    typedef decltype(T1()+(T2()*T3())+T2()) TRes;
//    Matrix<TRes> mt = mi + md + ms;
    
#ifndef NDEBUG
    std::cout << "Matrix<" << typeid(T1).name() << "> m1 = \n" << m1 << std::endl;
    std::cout << "Matrix<" << typeid(T2).name() << "> m2 = \n" << m2 << std::endl;
    std::cout << "Matrix<" << typeid(T3).name() << "> m3 = \n" << m3 << std::endl;
    std::cout << "Matrix<" << typeid(T4).name() << "> m4 = \n" << m4 << std::endl;

    std::cin.ignore();
#else
    std::cout << "m1 = \n" << m1 << std::endl;
    std::cout << "m2 = \n" << m2 << std::endl;
    std::cout << "m3 = \n" << m3 << std::endl;
    std::cout << "m4 = \n" << m4 << std::endl;
#endif

//    matrix::Matrix<decltype(T1()+T2()+T3())> mr = (m1+m2+m3);
    auto mr = m1+m2+m3+m4;
//    typedef decltype(m2*mr) MRes;

    CREATE_TEST(m2*mr)
    CREATE_TEST(m1+m2+m3+m4)
    CREATE_TEST((m1*m2)+(m3*m4))
    CREATE_TEST((m3+m4)*(m1+m2))
    CREATE_TEST(m1+m2*m3+m4)
    CREATE_TEST(m1+(m2+m3)+m4)
    CREATE_TEST(m1+(m2=(m3+m4)))
    CREATE_TEST(m2)

    return 0;
}
