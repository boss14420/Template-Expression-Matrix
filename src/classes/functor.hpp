/*
 * =====================================================================================
 *
 *       Filename:  functor.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/05/2011 06:18:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _FUNCTOR_HPP_
#define _FUNCTOR_HPP_

namespace boss14420 {
    namespace matrix {



#ifndef NDEBUG 
    #define MAKE_BINARY_OP_STRUCT(NAME, OP)  template<class T1, class T2> struct NAME { \
                typedef decltype(T1() OP T2()) T; \
                STRONG_INLINE T operator()(const T1& t1, const T2& t2) const { \
                    return t1 OP t2; \
                } \
                static std::string name() { \
                    return std::string( #NAME ) + "<" + typeid(T1).name() \
                        + ", " + typeid(T2).name() + ">"; \
                } \
            };
#else
    #define MAKE_BINARY_OP_STRUCT(NAME, OP)  template<class T1, class T2> struct NAME { \
                typedef decltype(T1() OP T2()) T; \
                STRONG_INLINE T operator()(const T1& t1, const T2& t2) const { \
                    return t1 OP t2; \
                } \
            };
#endif


            MAKE_BINARY_OP_STRUCT(scalar_sum_op, +)
            MAKE_BINARY_OP_STRUCT(scalar_diff_op, -)
            MAKE_BINARY_OP_STRUCT(scalar_xor_op, ^)
            MAKE_BINARY_OP_STRUCT(scalar_and_op, &)
            MAKE_BINARY_OP_STRUCT(scalar_or_op, |)


/*         template<class T1, class T2>
 *             struct scalar_sum_op {
 *                 typedef decltype(T1()+T2()) T;
 *                 EMPTY_STRUCT_CTOR(scalar_sum_op)
 * 
 *                 STRONG_INLINE T operator()(const T1& t1, const T2& t2) const {
 *                     return t1+t2;
 *                 }
 * 
 * #ifndef NDEBUG
 *                 static std::string name() {
 *                     return std::string("scalar_sum_op<") + typeid(T1).name()
 *                         + ", " + typeid(T2).name() + ">";
 * 
 *                 }
 * #endif
 *             };
 */



    }
}

#endif
