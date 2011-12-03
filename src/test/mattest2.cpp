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
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */


#include <sys/time.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iostream>

void test1(size_t n) {
    int* raw1 = new int[n*n];
    int* raw2 = new int[n*n];

    std::srand(std::time(NULL));
    for(size_t i = 0; i < n*n; ++i) {
        raw1[i] = std::rand() % 100;
        raw2[i] = std::rand() % 100;
    }
    
    int* rawr = new int[n*n];

    timeval begin, end;
    gettimeofday(&begin, NULL);

    int tmp;
    size_t idx1, idx2;
    for(size_t i=0; i<n; ++i) {
        idx1 = i * n;
        for(size_t j=0; j<n; ++j) {
            idx2 = j * n;
            tmp = 0;
            for(size_t k=0; k<n; ++k) {
                tmp += raw1[idx1+k] * raw2[idx2+k];
            }
            rawr[idx1+j] = tmp;
        }
    }

    gettimeofday(&end, NULL);
    double time = end.tv_sec - begin.tv_sec + 
        (end.tv_usec - begin.tv_usec) / 1000000.0;

    std::cout << "Test1, n = " << n << " : " << time << "s\n";

    delete[] raw1;
    delete[] raw2;
    delete[] rawr;

}

void test2(size_t n) {
    
    typedef int T1;
    typedef int T2;

    typedef int TMul;

    std::vector<int> raw1(n*n), raw2(n*n), rawr(n*n);

    std::srand(std::time(NULL));
    for(size_t i = 0; i < n*n; ++i) {
        raw1[i] = std::rand() % 100;
        raw2[i] = std::rand() % 100;
    }
    
    timeval begin, end;
    gettimeofday(&begin, NULL);

    TMul tmp;
    size_t idx1, idx2;
    for(size_t i=0; i<n; ++i) {
        idx1 = i * n;
        for(size_t j=0; j<n; ++j) {
            idx2 = j * n;
            tmp = 0;
            for(size_t k=0; k<n; ++k) {
                tmp += raw1[idx1+k] * raw2[idx2+k];
            }
            rawr[idx1+j] = tmp;
        }
    }

    gettimeofday(&end, NULL);
    double time = end.tv_sec - begin.tv_sec + 
        (end.tv_usec - begin.tv_usec) / 1000000.0;

    std::cout << "Test2, n = " << n << " : " << time << "s\n";

}

int main() {
    test1(1024);
    test2(1024);

    test1(2048);
    test2(2048);

    return 0;
}
