/*
 * =====================================================================================
 *
 *       Filename:  MatMul.hpp
 *
 *    Description:  Matrix multiplication
 *
 *        Version:  1.0
 *        Created:  11/28/2011 07:49:43 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef _MATMUL_HPP_
#define _MATMUL_HPP_

#include "matrix.hpp"
#ifndef NDEBUG_MT_PRODUCT_TIME
    #include <sys/time.h>
#endif


namespace boss14420 {
    namespace matrix {

        /*
         * =====================================================================================
         *        Class:  MatMul
         *  Description:  Matrix Multiplication Expression class
         * =====================================================================================
         */
        template < class T1, class T2>
            class MatMul : public MatExpression<decltype(T1()*T2()),RowMajor,MatMul<T1,T2>>
        {
            public:
                /* ====================  LIFECYCLE     ======================================= */
//                MatMul(const MatExpression<T1,E1>& u, const MatExpression<T2,E2>& v)    /* constructor */
//                    : _u(const_cast<MatExpression<T1,E1>&>(u).eval())
//                      , _v(const_cast<MatExpression<T2,E2>&>(v).eval())
//                template<MajorOrder mj1, MajorOrder mj2>
//                MatMul(Matrix<T1,mj1>& u, Matrix<T2,mj2>& v)    /* constructor */
                MatMul(const Matrix<T1,RowMajor>& u, const Matrix<T2, ColumnMajor>& v)
                    : _u(u), _v(v), _tmpres(u.row(), v.col()) {              

#ifndef NDEBUG_MT_PRODUCT_TIME
                        timeval begin, end;
                        gettimeofday(&begin, NULL);
#endif

                        const T1* raw1 = _u.getRawData();
                        const T2* raw2 = _v.getRawData();
                        TMul* rawr = _tmpres.getRawData();

                        TMul tmp;
                        size_t idx1, idx2;
                        for(size_t i=0; i<_u.row(); ++i) {
                            idx1 = i * _u.col();
                            for(size_t j=0; j<_v.col(); ++j) {
                                idx2 = j * _v.row();
                                tmp = 0;
                                for(size_t k=0; k<_u.col(); ++k) {
                                    tmp += raw1[idx1+k] * raw2[idx2+k];
                                }
                                rawr[idx1+j] = tmp;
                            }
                        }

#ifndef NDEBUG_MT_PRODUCT_TIME
                        gettimeofday(&end, NULL);
                        double time = end.tv_sec - begin.tv_sec + 
                            (end.tv_usec - begin.tv_usec) / 1000000.0;

                        std::cout << "Nhan ma tran " << _u.row() << "x" <<
                            _u.col() << " va " << _v.row() << "x" <<
                            _v.col() << " het : " << time << "s\n";
#endif

#ifndef NDEBUG
                          std::cout << "Constructor MatMul( "
//                              << MatExpression<T1,mj1,E1>::name() << ", "
//                              << MatExpression<T2,mj2,E2>::name() << " )\n";
                              << Matrix<T1,RowMajor>::name() << ", "
                              << Matrix<T2,ColumnMajor>::name() << ", )\n";
                          std::cout << "&u = " << &u << ", &v = " << &v << "\n";
                          std::cout << "&_u = " << &_u << ", &_v = " << &_v << "\n\n";
#endif

                      }

                MatMul(const MatMul<T1,T2>& mm) : _u(mm._u), _v(mm._v), _tmpres(mm._tmpres) {
#ifndef NDEBUG
                          std::cout << "Copy Constructor MatMul( " <<
                              name() << " )\n\n";
#endif
                }

                /* ====================  ACCESSORS     ======================================= */

                typedef decltype(T1() * T2()) TMul;

                size_t row() const { return _u.row();}
                size_t col() const { return _v.col();}
                size_t size() const { return _u.row() * _v.col();}

                TMul at(size_t i) const {
                    return _tmpres.at(i);
                }

                TMul at(size_t i, size_t j) const {
                    return _tmpres.at(i,j);
                }

                const Matrix<TMul,RowMajor>& eval() const {
                    return _tmpres;
                }

                /* ====================  MUTATORS      ======================================= */

                /* ====================  OPERATORS     ======================================= */

                TMul operator[](size_t i) const {
                    return _tmpres[i];
                }

                TMul operator()(size_t i, size_t j) const {
                    return _tmpres(i,j);
                }

                //        friend std::ostream& operator<< (std::ostream& os, MatAdd<T1,T2,E1,E2>& ma) {
                //            return os << ma.eval();
                //        }

#ifndef NDEBUG
                static std::string name() {
                    return std::string("MatMul<") + typeid(T1).name() + ", " +
                        typeid(T2).name() + ">";
                }
#endif

            protected:
                /* ====================  DATA MEMBERS  ======================================= */

            private:
                /* ====================  DATA MEMBERS  ======================================= */
                Matrix<T1,RowMajor> const& _u;
                Matrix<T2,ColumnMajor> const& _v;
                Matrix<TMul,RowMajor> _tmpres;

        }; /* ----------  end of template class MatMul ---------- */


        template <class T1, class T2, MajorOrder mj1, MajorOrder mj2, class E1, class E2>
            MatMul<T1, T2> // const
            operator* (MatExpression<T1,mj1,E1>& u, MatExpression<T2,mj2,E2>& v) {
                if(u.row() != v.col()) {
                    // throw
                }
                return MatMul<T1,T2>(u.eval(),v.eval());
            }
    }
}

#endif
