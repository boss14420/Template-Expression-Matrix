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

        /*
         * =====================================================================================
         *        Class:  MatPlus
         *  Description:  Matrix Add Expression class
         * =====================================================================================
         */
        template < class T1, class T2, MajorOrder mj1,class E1, class E2 >
            class MatPlus : public MatExpression<decltype(T1()+T2()),mj1,MatPlus<T1,T2,mj1,E1,E2> >
                           //class MatPlus : public MatExpression<decltype((T1)d1+(T2)d2), MatPlus<T1,T2,E1,E2> >
        {
            public:
                /* ====================  LIFECYCLE     ======================================= */
                MatPlus (MatExpression<T1,mj1,E1> const& u, MatExpression<T2,mj1,E2> const& v) 
                    : _u(u), _v(v), _tmpres(), _evaluted(false) {              /* constructor */
#ifndef NDEBUG
                          std::cout << "Constructor MatPlus ( "
                              << MatExpression<T1,mj1,E1>::name() << ", " 
                              << MatExpression<T2,mj1,E2>::name() << " )\n\n";
#endif

                      }

                MatPlus(const MatPlus<T1,T2,mj1,E1,E2>& mp) : _u(mp._u), _v(mp._v), _tmpres(mp._tmpres) {
#ifndef NDEBUG
                          std::cout << "Copy Constructor MatPlus ( "
                              << name() << " )\n\n";
#endif
                    }

                /* ====================  ACCESSORS     ======================================= */

                typedef decltype(T1() + T2()) TPlus;
                //        typedef decltype(T1(d1) + T2(d2)) TPlus;

                size_t row() const { return _u.row();}
                size_t col() const { return _u.col();}
                size_t size() const { return _u.size();}

                inline TPlus at(size_t i) const {
                    return _u.at(i) + _v.at(i);
                }

                inline TPlus at(size_t i, size_t j) const {
                    return _u.at(i,j) + _v.at(i,j);
                }

                const Matrix<TPlus,mj1>& eval() {
                    if(!_evaluted) {
                        for(size_t i = 0; i < size(); ++i)
                            _tmpres[i] = _u[i] + _v[i];

                    }
                    return _tmpres;
                }

#ifndef NDEBUG
                static std::string name() {
                    return std::string("MatPlus<") + typeid(T1).name() + ", " +
                        typeid(T2).name() + ", " 
                        + ((mj1 == RowMajor) ? "RowMajor, " : "ColumnMajor, ")
                        + E1::name() + ", " + E2::name() + ">";
                }
#endif

                /* ====================  MUTATORS      ======================================= */

                /* ====================  OPERATORS     ======================================= */

                inline TPlus operator[](size_t i) const {
                    return _u[i] + _v[i];
                }

                inline TPlus operator()(size_t i, size_t j) const {
                    return _u(i,j) + _v(i,j);
                }

                //        MatExpression<TPlus, MatPlus<T2,T2,E1,E2>>& operator MatExpression<TPlus, MatPlus<T2,T2,E1,E2>>& () {
                //            return static_cast<MatExpression<TPlus, MatPlus<T2,T2,E1,E2>>&>(*this);
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
                Matrix<TPlus,mj1> _tmpres;
                bool _evaluted;

        }; /* ----------  end of template class MatPlus  ---------- */



        template <class T1, class T2, MajorOrder mj1, class E1, class E2>
            inline MatPlus<T1, T2, mj1, E1, E2> //const
            operator+ (MatExpression<T1,mj1,E1> const& u, MatExpression<T2,mj1,E2> const& v) {
                if(u.row() != v.row() || u.col() != v.col()) {
                    // throw 
                }
                return MatPlus<T1,T2,mj1,E1,E2>(u,v);
            }

        template <class T1, class T2, class E1, class E2>
            inline MatPlus<T1, T2, RowMajor, E1, Matrix<T2,RowMajor>>
            operator+ (MatExpression<T1,RowMajor,E1> const& u, MatExpression<T2,ColumnMajor,E2>& v) {
                if(u.row() != v.row() || u.col() != v.col()) {
                    // throw
                }
                Matrix<T2,RowMajor> vv(v.eval());
                return MatPlus<T1,T2,RowMajor,E1,Matrix<T2,RowMajor>>(u,vv);
            }

        template <class T1, class T2, class E1, class E2>
            MatPlus<T1, T2, ColumnMajor, E1, Matrix<T2,ColumnMajor>>
            inline operator+ (MatExpression<T1,ColumnMajor,E1> const& u, MatExpression<T2,RowMajor,E2>& v) {
                if(u.row() != v.row() || u.col() != v.col()) {
                    // throw
                }
                Matrix<T2,ColumnMajor> vv(v.eval());
                return MatPlus<T1,T2,ColumnMajor,E1,Matrix<T2,ColumnMajor>>(u,vv);
            }

    }
}

#endif
