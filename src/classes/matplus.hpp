/*
 * =====================================================================================
 *
 *       Filename:  MatPlus.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/28/2011 07:53:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _MATPLUS_HPP_
#define _MATPLUS_HPP_

#include "matrix.hpp"

namespace boss14420 {
    namespace matrix {

        template<class _E1, class _E2>
            struct traits<MatPlus<_E1,_E2>> {
                typedef _E1 E1;
                typedef _E2 E2;

                typedef typename traits<E1>::T T1;
                typedef typename traits<E2>::T T2;

                typedef decltype(T1()+T2()) T;
                static const MajorOrder major_order = traits<E1>::major_order;
            };

        /*
         * =====================================================================================
         *        Class:  MatPlus
         *  Description:  Matrix Add Expression class
         * =====================================================================================
         */
        template < class E1, class E2 >
            class MatPlus : public MatExpression<MatPlus<E1,E2>>
//            class MatPlus : public MatExpression<decltype(T1()+T2()),major_order,MatPlus<T1,T2,major_order,E1,E2>>
        {
            public:
                typedef typename traits<E1>::T T1;
                typedef typename traits<E2>::T T2;

                static const MajorOrder major_order = traits<E1>::major_order;
                typedef decltype(T1() + T2()) T;


                /* ====================  LIFECYCLE     ======================================= */
                MatPlus (MatExpression<E1> const& u, MatExpression<E2> const& v) 
                    : _u(u), _v(v), _tmpres(), _evaluted(false) {              /* constructor */
#ifndef NDEBUG
                          std::cout << "Constructor MatPlus ( "
                              << MatExpression<E1>::name() << ", " 
                              << MatExpression<E2>::name() << " )\n\n";
#endif

                      }

                MatPlus(const MatPlus<E1,E2>& mp) : _u(mp._u), _v(mp._v), _tmpres(mp._tmpres) {
#ifndef NDEBUG
                          std::cout << "Copy Constructor MatPlus ( "
                              << name() << " )\n\n";
#endif
                    }

                /* ====================  ACCESSORS     ======================================= */

                size_t row() const { return _u.row();}
                size_t col() const { return _u.col();}
                size_t size() const { return _u.size();}

                STRONG_INLINE T at(size_t i) const {
                    return _u.at(i) + _v.at(i);
                }

                STRONG_INLINE T at(size_t i, size_t j) const {
                    return _u.at(i,j) + _v.at(i,j);
                }

                const Matrix<T,major_order>& eval() {
                    if(!_evaluted) {
                        size_t i;
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
                        for(i = 0; i < size(); ++i)
                            _tmpres[i] = _u[i] + _v[i];

                    }
                    return _tmpres;
                }

#ifndef NDEBUG
                static std::string name() {
//                    return std::string("MatPlus<") + typeid(T1).name() + ", " +
//                        typeid(T2).name() + ", " 
//                        + ((major_order == RowMajor) ? "RowMajor, " : "ColumnMajor, ")
//                        + E1::name() + ", " + E2::name() + ">";
                    return std::string("MatPlus<") + E1::name() + E2::name() + ">";
                }
#endif

                /* ====================  MUTATORS      ======================================= */

                /* ====================  OPERATORS     ======================================= */

                STRONG_INLINE T operator[](size_t i) const {
                    return _u[i] + _v[i];
                }

                STRONG_INLINE T operator()(size_t i, size_t j) const {
                    return _u(i,j) + _v(i,j);
                }

                //        MatExpression<T, MatPlus<T2,T2,E1,E2>>& operator MatExpression<T, MatPlus<T2,T2,E1,E2>>& () {
                //            return static_cast<MatExpression<T, MatPlus<T2,T2,E1,E2>>&>(*this);
                //        }

                //        friend std::ostream& operator<< (std::ostream& os, MatPlus<T1,T2,E1,E2>& ma) {
                //            return os << ma.eval();
                //        }

            protected:
                /* ====================  DATA MEMBERS  ======================================= */

            private:
                /* ====================  DATA MEMBERS  ======================================= */
                E1 const& _u;
                E2 const& _v;
                Matrix<T,major_order> _tmpres;
                bool _evaluted;

        }; /* ----------  end of template class MatPlus  ---------- */



        template <class E1, class E2>
            STRONG_INLINE MatPlus<E1, E2> //const
            operator+ (MatExpression<E1> const& u, MatExpression<E2> const& v) {
                if(u.row() != v.row() || u.col() != v.col()) {
                    // throw 
                }

                return MatPlus<E1,E2>(u,v);
            }

    }
}

#endif
