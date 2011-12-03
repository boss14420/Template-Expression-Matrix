/*
 * =====================================================================================
 *
 *       Filename:  random.hpp
 *
 *    Description:  use to generate random value
 *
 *        Version:  1.0
 *        Created:  11/28/2011 12:00:18 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <cstdlib>
#include <limits>
#include <ctime>

namespace boss14420 {
    namespace random {

        struct seed {
            seed() {
                srand(std::time(NULL));
            }
        };

        seed sd;

        template<class T, bool specialized = std::numeric_limits<T>::is_specialized>
            class random {
                const T& _min, _max;
                T _gap;
              public:
                random( const T& min=std::numeric_limits<T>::min(),
                        const T& max=std::numeric_limits<T>::max()) 
                            : _min(min), _max(max) {
                                
                                if(max - min < RAND_MAX)
                                    _gap = max - min +1;
                                else
                                    _gap = max - min;
                }

                T nextRan() {
                    if(_gap < RAND_MAX) {
                        unsigned tmp = std::rand();
                        T mod = tmp - unsigned(tmp / (_gap*5)) * (_gap*5);
                        return mod / 5 + _min;
                    } else {
                        return std::rand()*(_gap/RAND_MAX) + std::rand() + _min;
                    }
                }
            };

        template<class T>
            class random<T, false> {
                const T& _min, _max;

              public:
                random( const T& min, const T& max) : _min(min), _max(max) {
                }

                T nextRan() {
                    return T::random(_min, _max);
                }
            };

    }
}

#endif
