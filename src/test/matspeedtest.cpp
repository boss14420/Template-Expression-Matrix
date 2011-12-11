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
#include <sys/time.h>

#ifndef NDEBUG
    #define NDEBUG 1
#endif

#include "../classes/matrix"
#include "../trivial/matrix.hpp"
#include <eigen3/Eigen/Dense>


void test_product_time(size_t n, size_t loops) {
    
    typedef float T;

    using namespace boss14420;

#ifndef NDEBUG_MT_PRODUCT_TIME
//    typedef decltype(T1()*T2()) TRes;

    timeval begin, end;
    double time;

    // boss14420 matrix
    matrix::Matrix<T> m1 = matrix::Matrix<T>::Random(n,n,-100,100);
    matrix::Matrix<T,matrix::ColumnMajor> m2 = matrix::Matrix<T,matrix::ColumnMajor>::Random(n,n,-100,100);

    matrix::Matrix<T> m3;
    gettimeofday(&begin, NULL);
    for(size_t i = 0; i < loops; ++i)
        m3 = m1 * m2;
    gettimeofday(&end, NULL);

    time = end.tv_sec - begin.tv_sec + 
        (end.tv_usec - begin.tv_usec) / 1000000.0;

    std::cout << "boss14420 matrix: Thoi gian nhan 2 ma tran " << n << "x" << n 
        << ", " << loops << " lan la : " << time << std::endl;
//
//    // Eigei matrix
//
//    using Eigen::Dynamic;
//
//    typedef Eigen::Matrix<T,Dynamic,Dynamic> EigenMatDD;
//
//    EigenMatDD em1 = EigenMatDD::Random(n,n);
//    EigenMatDD em2 = EigenMatDD::Random(n,n);
//
//    EigenMatDD emr;
//    gettimeofday(&begin, NULL);
//    for(size_t i = 0; i < loops; ++i)
//        emr = em1*em2;
//    gettimeofday(&end, NULL);
//
//    time = end.tv_sec - begin.tv_sec + 
//        (end.tv_usec - begin.tv_usec) / 1000000.0;
//
//    std::cout << "Eigen matrix: Thoi gian nhan 2 ma tran " << n << "x" << n 
//        << ", " << loops << " lan la : " << time << std::endl << std::endl;
//
#endif
}

void test_expression_time(size_t n, size_t loops) {

#ifndef NDEBUG_MT_EXPRESSION_TIME
    typedef float T;
    const size_t SZ = 9;

    using namespace boss14420;
    using Eigen::Dynamic;

    matrix::Matrix<T> m[SZ];
//    trivial::Matrix<T> tm[SZ];
    Eigen::Matrix<T,Dynamic,Dynamic> em[SZ];

    for(size_t i = 0; i < SZ; ++i) {
        m[i] = matrix::Matrix<T>::Random(n,n,0,100);
//        tm[i] = trivial::Matrix<T>::Random(n,n,0,100);
        em[i] = Eigen::Matrix<T,Dynamic,Dynamic>::Random(n,n);
    }

    timeval begin, end;
    double time;

    // Optimized matrix
//    matrix::Matrix<T> *mr = new matrix::Matrix<T>[loops];

    gettimeofday(&begin, NULL);

    for(size_t i = 0; i < loops; ++i) {
//        mr[i] = m[0] + m[1] + m[2] + m[3] + m[4]
//            + m[5] + m[6] + m[7] + m[8] + m[9];
    matrix::Matrix<T> mr = m[std::rand()%SZ] 
//        val = (m[std::rand()%SZ] + m[std::rand()%SZ] 
            +m[std::rand()%SZ] +m[std::rand()%SZ] +m[std::rand()%SZ]
            +m[std::rand()%SZ] +m[std::rand()%SZ] +m[std::rand()%SZ]
            +m[std::rand()%SZ] +m[std::rand()%SZ] +m[std::rand()%SZ];
    }

    gettimeofday(&end, NULL);

    time = end.tv_sec - begin.tv_sec + 
        (end.tv_usec - begin.tv_usec) / 1000000.0;

    std::cout << "template expression: Thoi gian cong " 
        << SZ << " ma tran " << n << "x" << n << ", " 
        << loops << " lan la : " << time << std::endl;

//    delete[] mr;
//    mr = NULL;


    // Eigen matrix

    gettimeofday(&begin, NULL);

    for(size_t i = 0; i < loops; ++i) {
        Eigen::Matrix<T,Eigen::Dynamic,Dynamic> emr = em[std::rand()%SZ] 
            +em[std::rand()%SZ] +em[std::rand()%SZ] +em[std::rand()%SZ]
            +em[std::rand()%SZ] +em[std::rand()%SZ] +em[std::rand()%SZ]
            +em[std::rand()%SZ] +em[std::rand()%SZ] +em[std::rand()%SZ];
    }

    gettimeofday(&end, NULL);

    time = end.tv_sec - begin.tv_sec + 
        (end.tv_usec - begin.tv_usec) / 1000000.0;

    std::cout << "Eigen Matrix : Thoi gian cong " 
        << SZ << " ma tran " << n << "x" << n << ", " 
        << loops << " lan la : " << time << std::endl;

//    delete[] mr;
//    mr = NULL;




    // trivial matrix
//    trivial::Matrix<T> *tmr = new trivial::Matrix<T>[loops];
    
//    gettimeofday(&begin, NULL);
//
//    for(size_t i = 0; i < loops; ++i) {
////        tmr[i] = tm[0] + tm[1] + tm[2] + tm[3] + tm[4]
////            + tm[5] + tm[6] + tm[7] + tm[8] + tm[9];
//    trivial::Matrix<T> tmr = tm[std::rand()%SZ] 
////        val = (tm[std::rand()%SZ] + tm[std::rand()%SZ] 
//            +tm[std::rand()%SZ] +tm[std::rand()%SZ] +tm[std::rand()%SZ]
//            +tm[std::rand()%SZ] +tm[std::rand()%SZ] +tm[std::rand()%SZ]
//            +tm[std::rand()%SZ] +tm[std::rand()%SZ] +tm[std::rand()%SZ];
//    }
//
//    gettimeofday(&end, NULL);
//
//    time = end.tv_sec - begin.tv_sec + 
//        (end.tv_usec - begin.tv_usec) / 1000000.0;
//
//    std::cout << "trivial matrix: Thoi gian cong " 
//        << SZ << " ma tran " << n << "x" << n << ", " 
//        << loops << " lan la : " << time << std::endl << std::endl;

//    delete[] tmr;
//    tmr = NULL;

#endif

}

int main(int argc, char *argv[]) {

    char mode = '0';
    if(argc >= 2)
        mode = argv[1][0];

    if('0' == mode) {
        test_product_time(100,4000);
        test_product_time(200,1000);
        test_product_time(512,200);
        test_product_time(1024,1);
        test_product_time(2048,1);
        test_product_time(4096,1);
    }

    if('1' == mode) {
        test_expression_time(4,800000);
        test_expression_time(8,400000);
        test_expression_time(16,200000);
        test_expression_time(64,100000);
        test_expression_time(128,50000);
        test_expression_time(256,2000);
        test_expression_time(512,4000);
        test_expression_time(1024,1000);
        test_expression_time(4096,500);
    }

    return 0;
}

