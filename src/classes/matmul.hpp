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

#ifdef __INTEL_COMPILER
    #include <boost/shared_ptr.hpp>
    using boost::shared_ptr;
#else
    #include <memory>
    using std::shared_ptr;
#endif

#include <utility>
#include "matrix.hpp"
//#ifndef NDEBUG_MT_PRODUCT_TIME
//    #include <sys/time.h>
//#endif


namespace boss14420 {
    namespace matrix {

        template<class _E1, class _E2>
            struct traits<MatMul<_E1,_E2>> {
                typedef _E1 E1;
                typedef _E2 E2;

                typedef typename traits<E1>::T T1;
                typedef typename traits<E2>::T T2;
                typedef decltype(T1()*T2()+T1()*T2()) T;

//                typedef Matrix<T1,RowMajor> ValueE1;
//                typedef Matrix<T2,ColumnMajor> ValueE2;
                typedef typename nested<E1,NeedEvalute |
                    ((traits<E1>::major_order == ColumnMajor) ? NeedTranposeMask : 0) >::type ValueE1;
                typedef typename nested<E2,NeedEvalute |
                    ((traits<E2>::major_order == RowMajor) ? NeedTranposeMask : 0)>::type ValueE2;

                typedef typename nested<E1,0>::type NestedE1;
                typedef typename nested<E2,0>::type NestedE2;

                typedef Matrix<T,RowMajor> ReturnType;

                static const MajorOrder major_order = RowMajor;
            };

        template<class E1, class E2> 
            struct nested<MatMul<E1,E2>,NeedTranposeMask|NeedEvalute> {
                typedef Matrix<typename traits<MatMul<E1,E2>>::T,ColumnMajor> type;
        };

        template<class E1, class E2, int Option>
            struct nested<MatMul<E1,E2>,Option> {
//                typedef Matrix<typename traits<MatMul<E1,E2>>::T, RowMajor> type;
                typedef MatMul<E1,E2> type;
            };

//        template<class T>
//            struct pass_to_expression<Matrix<T, RowMajor>> {
//                template<class E1, class E2>
//                Matrix<T,RowMajor> 
//                    operator()(const MatMul<E1,E2>& mm) {
//                        return std::move((const_cast<MatMul<E1,E2>&>(mm))._tmpres);
//                    }
//
//                template<class E1, class E2>
//                Matrix<T,RowMajor>
//                    operator()(const MatExpression<MatMul<E1,E2>>& mm) {
//                        return std::move(((MatMul<E1,E2>)(const_cast<MatExpression<MatMul<E1,E2>>&>(mm)))._tmpres);
//                    }
//            };

        /*
         * =====================================================================================
         *        Class:  MatMul
         *  Description:  Matrix Multiplication Expression class
         * =====================================================================================
         */
        template < class E1, class E2>
            class MatMul : public MatExpression<MatMul<E1,E2>>
        {
            public:

                typedef typename traits<E1>::T T1;
                typedef typename traits<E2>::T T2;
                typedef decltype(T1()*T2()+T1()*T2()) T;

//                typedef Matrix<T1,RowMajor> ValueE1;
//                typedef Matrix<T2,ColumnMajor> ValueE2;
                typedef typename nested<E1,NeedEvalute |
                    ((traits<E1>::major_order == ColumnMajor) ? NeedTranposeMask : 0) >::type ValueE1;
                typedef typename nested<E2,NeedEvalute |
                    ((traits<E2>::major_order == RowMajor) ? NeedTranposeMask : 0)>::type ValueE2;

                typedef typename nested<E1,0>::type NestedE1;
                typedef typename nested<E2,0>::type NestedE2;

                typedef Matrix<T,RowMajor> ReturnType;

                static const MajorOrder major_order = RowMajor;

//                friend template class ReturnType;
//                friend struct pass_to_expression<ReturnType>;


                /* ====================  LIFECYCLE     ======================================= */
//                MatMul(const MatExpression<T1,E1>& u, const MatExpression<T2,E2>& v)    /* constructor */
//                    : _u(const_cast<MatExpression<T1,E1>&>(u).eval())
//                      , _v(const_cast<MatExpression<T2,E2>&>(v).eval())
//                template<MajorOrder mj1, MajorOrder mj2>
//                MatMul(Matrix<T1,mj1>& u, Matrix<T2,mj2>& v)    /* constructor */
                MatMul(const MatExpression<E1>& u, const MatExpression<E2>& v)
                    : _u(pass_to_expression<NestedE1>()(u)),
                      _v(pass_to_expression<NestedE2>()(v)),
//                    : _u(u),
//                      _v(v),
                      _evaluted(false), _tmpres(u.row(),v.col()) {

                        _evaluted = true;
//                        eval();
                        
                        //#ifndef NDEBUG_MT_PRODUCT_TIME
                        //timeval begin, end;
                        //gettimeofday(&begin, NULL);
                        //#endif

                        const size_t m = _u.row(), n = _u.col(), p = _v.col();
//                        Matrix<T,RowMajor> _tmpres(m,p);

                        ValueE1 vu = ValueE1(_u);
                        ValueE2 vv = ValueE2(_v);

                        T1 const *raw1 = vu.getRawData();
                        T2 const *raw2 = vv.getRawData();
                        T* rawr = _tmpres.getRawData();

                        T tmp;
                        size_t i,j,k,idx1=0,idx2,idx3=0;
                        //const T1* rIndex1 = raw1;
                        //const T2* rIndex2;
                        //T* rIndexr = rawr;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k,tmp,idx1,idx2,idx3)
#endif
                        for(i=0; i<m; ++i) {
                            idx1 = i * n;
                            idx3 = i * p;
                            idx2 = 0;
                            for(j=0; j<p; ++j) {
                                tmp = 0;
                                for(k=0; k<n; ++k) {
                                    tmp += raw1[idx1+k] * raw2[idx2+k];
                                }
                                rawr[idx3+j] = tmp;
                                idx2 += n;
                            }
                            //                            idx1 += n;
                            //                            idx3 += p;
                        }

                        raw1 = NULL;
                        raw2 = NULL;
                        rawr = NULL;
                        

                        /*                         for(i=0; i<m; ++i) {
                         *                             rIndex2 = raw2;
                         *                             for(j=0; j<p; ++j) {
                         *                                 tmp = 0;
                         *                                 for(k=0; k<n; ++k) {
                         *                                     tmp += rIndex1[k] * rIndex2[k];
                         *                                 }
                         *                                 rIndexr[j] = tmp;
                         *                                 rIndex2 += n;
                         *                             }
                         *                             rIndex1 += n;
                         *                             rIndexr += p;
                         *                         }
                         */

                        //#ifndef NDEBUG_MT_PRODUCT_TIME
                        //                        gettimeofday(&end, NULL);
                        //                        double time = end.tv_sec - begin.tv_sec + 
                        //                            (end.tv_usec - begin.tv_usec) / 1000000.0;
                        //
                        //                        std::cout << "Nhan ma tran " << _u.row() << "x" <<
                        //                            _u.col() << " va " << _v.row() << "x" <<
                        //                            _v.col() << " het : " << time << "s\n";
                        //#endif


#ifndef NDEBUG
                          std::cout << "Constructor MatMul( "
//                              << MatExpression<T1,mj1,E1>::name() << ", "
//                              << MatExpression<T2,mj2,E2>::name() << " )\n";
                              << MatExpression<E1>::name() << ", "
                              << MatExpression<E2>::name() << " ) at " << this << " \n";
                          std::cout << "&u = " << &u << ", &v = " << &v << "\n";
//                          std::cout << "u = " << u << "\nv = " << v << "\n";
                          std::cout << "&_tmpres = " << &_tmpres << ", _tmpres = " << _tmpres << "\n\n";
#endif

                      }

                MatMul(MatMul<E1,E2>&& mm) : 
                    _u(pass_to_expression<NestedE1>()(mm._u)), 
                    _v(pass_to_expression<NestedE2>()(mm._v)), 
                    _evaluted(mm._evaluted), 
                    _tmpres(std::move(mm._tmpres)) {
#ifndef NDEBUG
                          std::cout << "Move Constructor MatMul( " <<
                              name() << " )\n\n";
#endif
                }

                /* ====================  ACCESSORS     ======================================= */

                STRONG_INLINE size_t row() const { return _tmpres.row();}
                STRONG_INLINE size_t col() const { return _tmpres.col();}
                STRONG_INLINE size_t size() const { return _tmpres.size();}

                STRONG_INLINE T at(size_t i) const {
                    return _tmpres.at(i);
                }

                STRONG_INLINE T at(size_t i, size_t j) const {
                    return _tmpres.at(i,j);
                }

                void eval() {

//                    if(!_evaluted) {
//
//                                                }
//
//                        _evaluted = true;
//                    }

                }

                /* ====================  MUTATORS      ======================================= */

                /* ====================  OPERATORS     ======================================= */

                STRONG_INLINE T operator[](size_t i) const {
                    return _tmpres[i];
                }

                STRONG_INLINE T operator()(size_t i, size_t j) const {
                    return _tmpres(i,j);
                }

                //        friend std::ostream& operator<< (std::ostream& os, MatAdd<T1,T2,E1,E2>& ma) {
                //            return os << ma.eval();
                //        }

#ifndef NDEBUG
                static std::string name() {
                    return std::string("MatMul<") + E1::name() + ", " +
                        E2::name() + ">";
                }
#endif

                operator ReturnType () const {
                    return std::move(_tmpres);
                }

                operator const ReturnType& () const {
                    return std::move(_tmpres);
                }

            protected:
                /* ====================  DATA MEMBERS  ======================================= */

            private:
                /* ====================  DATA MEMBERS  ======================================= */
//                Matrix<T1,RowMajor> const& _u;
//                Matrix<T2,ColumnMajor> const& _v;
                NestedE1 _u;
                NestedE2 _v;

//                ReturnType _tmpres;
                bool _evaluted;
//                shared_ptr<ReturnType> _ptrres;
                ReturnType _tmpres;

        }; /* ----------  end of template class MatMul ---------- */


//        template <class E1, class E2>
//            MatMul<typename traits<E1>::T, typename traits<E2>::T> // const
//            operator* (MatExpression<E1> const& u, MatExpression<E2> const& v) {
//                if(u.row() != v.col()) {
//                    // throw
//                }
//                typedef typename traits<E1>::T T1;
//                typedef typename traits<E2>::T T2;
//
//                MatExpression<E1>& _u = const_cast<MatExpression<E1>&>(u);
//                MatExpression<E2>& _v = const_cast<MatExpression<E2>&>(v);
//
////                MatExpression<E1>& _u = u;
////                MatExpression<E2>& _v = v;
//
//                if(traits<E1>::major_order == RowMajor) {
//                    if(traits<E2>::major_order == RowMajor)
//                        return MatMul<T1,T2>(_u.eval(), Matrix<T2,ColumnMajor>(_v));
//                    else
//                        return MatMul<T1,T2>(_u.eval(), _v.eval());
//                } else {
//                    if(traits<E2>::major_order == RowMajor)
//                        return MatMul<T1,T2>(Matrix<T1,RowMajor>(_u),Matrix<T2,ColumnMajor>(_v));
//                    else
//                        return MatMul<T1,T2>(Matrix<T1,RowMajor>(_u),_v.eval());
//                }
//            }

        template <class E1, class E2>
            MatMul<E1,E2> operator* (MatExpression<E1> const& u, MatExpression<E2> const& v) {
                return MatMul<E1,E2>(u,v);
            }
    }
}

#endif
