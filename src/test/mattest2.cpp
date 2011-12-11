/*
 * =====================================================================================
 *
 *       Filename:  mattest2.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/02/2011 01:02:26 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS15520 (boss15520), boss15520@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */


#include <sys/time.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iostream>
#include "../classes/matrix"
#include "../trivial/matrix.hpp"
#include <eigen3/Eigen/Dense>

#ifndef __TEST_TYPE__
#define __TEST_TYPE__ float
#endif

int main(int argc, char *argv[]) {
    using namespace boss14420;

    size_t SIZE = 4096;

    if(argc == 2)
        SIZE = std::atoi(argv[1]);

    std::cout << "SIZE = " << SIZE << std::endl;

    timeval begin, end;
    double time;

    {
        // Trivial Matrix
        trivial::Matrix<__TEST_TYPE__> tm[5];

        std::cout << "Generate random matrix ..." << std::endl;
        for(int i = 0; i < 4; ++i) {
            tm[i] = trivial::Matrix<__TEST_TYPE__>::Random(SIZE,SIZE,0,100);
        }

        std::cout << "trivial matrix: " ;
        std::cout.flush();
        gettimeofday(&begin, NULL);
        tm[4] = tm[0] + tm[1] + tm[2] + tm[3];
        gettimeofday(&end, NULL);

        time = end.tv_sec - begin.tv_sec + 
            (end.tv_usec - begin.tv_usec) / 1000000.0;

        std::cout << time << std::endl;

    }

    {
        // Template Expression

        matrix::Matrix<__TEST_TYPE__> mm[5];
        std::cout << "Generate random matrix ..." << std::endl;
        for(int i = 0; i < 4; ++i) {
            mm[i] = matrix::Matrix<__TEST_TYPE__>::Random(SIZE,SIZE,0,100);
        }

        std::cout << "template expression matrix: ";
        std::cout.flush();
        gettimeofday(&begin, NULL);
        mm[4] = mm[0] + mm[1] + mm[2] + mm[3];
        gettimeofday(&end, NULL);

        time = end.tv_sec - begin.tv_sec + 
            (end.tv_usec - begin.tv_usec) / 1000000.0;

        std::cout << time << std::endl;
    }

    {
        // Eigen Matrix

        using Eigen::Dynamic;
        Eigen::Matrix<__TEST_TYPE__,Dynamic, Dynamic> em[5];

        std::cout << "Generate random matrix ..." << std::endl;
        for(int i = 0; i < 4; ++i) {
            em[i] = Eigen::Matrix<__TEST_TYPE__,Dynamic,Dynamic>::Random(SIZE,SIZE);
        }

        std::cout << "eigen matrix: ";
        std::cout.flush();
        gettimeofday(&begin, NULL);
        em[4] = em[0] + em[1] + em[2] + em[3];
        gettimeofday(&end, NULL);

        time = end.tv_sec - begin.tv_sec + 
            (end.tv_usec - begin.tv_usec) / 1000000.0;

        std::cout << time << std::endl;
    }

    return 0;
}
