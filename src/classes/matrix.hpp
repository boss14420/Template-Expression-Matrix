/*
 * =====================================================================================
 *
 *       Filename:  matrix.hpp
 *
 *    Description:  matrix class
 *
 *        Version:  1.0
 *        Created:  11/27/2011 06:59:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

//#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <type_traits>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <cstring>

#include "random.hpp"


// GCC has not implement this
template<class T>
struct is_trivially_copyable {
    static const bool value = true;
};

namespace boss14420 {
    namespace matrix {

//        enum MajorOrder { RowMajor = 0, ColumnMajor = 1 };
        typedef bool MajorOrder;
        const bool RowMajor = true;
        const bool ColumnMajor = false;

        const int NeedEvalute = 1;
        const int NeedTranposeMask = 2;
        
        template<class Derived> struct traits;

        template<class T, MajorOrder mj> class Matrix;
        template<class BinaryOp, class E1, class E2> class CWiseBinaryOp;
        template<class E1, class E2> class CWiseBinaryOp2;
        template<class E1, class E2> class MatPlus;
        template<class E1, class E2> class MatMul;

        template<class Derived, int Option> struct nested;

        /*
         * =====================================================================================
         *        Class:  MatExpression
         *  Description:  template Matrix expression class
         * =====================================================================================
         */
        template <class E>
            class MatExpression
            {
                public:

                    typedef typename std::remove_reference<E>::type noref_E;
                    typedef typename std::remove_const<noref_E>::type norefconst_E;
                    typedef typename traits<norefconst_E>::T T;
                    static const MajorOrder major_order = traits<norefconst_E>::major_order;
                    static const size_t ParallelThreshold = 4096 / sizeof(T);

                    /* ====================  LIFECYCLE     ======================================= */

                    /* ====================  ACCESSORS     ======================================= */
                    size_t row() const {
                        return ((E const&)*this).row();
                    }

                    size_t col() const {
                        return ((E const&)*this).col();
                    }

                    size_t size() const {
                        return ((E const&)*this).size();
                    }

                    T at(size_t i) const {
                        return ((E const&)*this).at(i);
                    }

                    T at(size_t i, size_t j) const {
                        return ((E const&)*this).at(i,j);
                    }

#ifndef NDEBUG
                    static std::string name() {
                        //                        return std::string("MatExpression<") + typeid(T).name() 
                        //                            + ", " + E::name() + "> ";
                        return E::name();
                    }
#endif

                    /* ====================  MUTATORS      ======================================= */
                    const Matrix<T, major_order>& eval() {
                        return static_cast<E&>(*this).eval();
                    }

                    /* ====================  OPERATORS     ======================================= */

                    STRONG_INLINE T operator[](size_t i) const {
                        return ((E const&)*this)[i];
                    }

                    STRONG_INLINE T operator()(size_t i, size_t j) const {
                        return ((E const&)*this)(i, j);
                    }

#ifdef __INTEL_COMPILER
                    STRONG_INLINE operator E&() {
                        return reinterpret_cast<E&>(*this);
                    }

                    STRONG_INLINE operator E const&() const {
                        return reinterpret_cast<E const&>(*this);
                    }
#else
                    STRONG_INLINE operator E&() {
                        return static_cast<E&>(*this);
                    }

                    STRONG_INLINE operator E const&() const {
                        return static_cast<E const&>(*this);
                    }
#endif
//                    friend std::ostream& operator<< (std::ostream& os, const MatExpression<E>& me) {
//                        //            if(typeid(E) != typeid(Matrix<T>)) {
//                        //                return os << static_cast<E&>(me).eval();
//                        //            } else {
//                        //                return os << static_cast<Matrix<T> const&>(me);
//                        //            }
//                        return os << static_cast<Matrix<T,major_order> const&>(me);
//                    }

                    friend std::ostream& operator<< (std::ostream& os, const MatExpression<E>& mt) {

//                        E const& mt = ((E const&)me);

                        if(&os == &std::cout || !is_trivially_copyable<T>::value) { // standard input
                            os << mt.row() << "x" << mt.col() << std::endl;

                            for(size_t i = 0; i < mt.row(); ++i) {
                                for(size_t j = 0; j < mt.col(); ++j) {
                                    os << mt(i,j) << " ";
                                }
                                os << std::endl;
                            }
                        } else {
                            char mj = (char)major_order;
                            size_t row = mt.row();
                            size_t col = mt.col();
                            os.write((char*)&row, sizeof(size_t))
                                .write((char*)&col, sizeof(size_t))
                                .write((char*)&mj, sizeof(MajorOrder));

                            T tmp;
                            for(size_t i = 0; i < mt.size(); ++i) {
                                tmp = mt[i];
                                os.write((char*)&tmp, sizeof(T));
                            }
                            //                std::ostreambuf_iterator<T> obi(os);
                            //                std::copy(mt._data.begin(), mt._data.end(), obi);
                        }

                        return os;
                    }

                protected:
                    /* ====================  DATA MEMBERS  ======================================= */

                private:
                    /* ====================  DATA MEMBERS  ======================================= */

            }; /* ----------  end of template class MatExpression  ---------- */


        template<class _T>
            struct traits<Matrix<_T,RowMajor>> {
                typedef _T T;
                static const MajorOrder major_order = RowMajor;
            };

        template<class _T>
            struct traits<Matrix<_T,ColumnMajor>> {
                typedef _T T;
                static const MajorOrder major_order = ColumnMajor;
            };

        template<class T, MajorOrder major_order>
            struct nested<Matrix<T,major_order>,3> {
                typedef Matrix<T, !major_order> type;
            };

        template<class T, MajorOrder major_order, int Option>
            struct nested<Matrix<T,major_order>,Option> {
                typedef const Matrix<T,major_order>& type;
            };

        template<class Derived>
            struct pass_to_expression {
                Derived&& operator()(const MatExpression<Derived>& d) {
                    return std::move((Derived&)(const_cast<MatExpression<Derived>&>(d)));
                }

                Derived&& operator()(const Derived& d) {
                    return std::move(const_cast<Derived&>(d));
                }
            };
    
//        template<class E>
//            struct pass_to_expression {
//                E operator()(MatExpression<E>& me) {
//                    return me;
//                }
//
//                E operator()(E& e) {
//                    return e;
//                }
//            };

        template<class T, MajorOrder major_order>
            struct pass_to_expression<const MatExpression<Matrix<T,major_order>>&> {
                const Matrix<T,major_order>& 
                    operator()(const MatExpression<Matrix<T,major_order>>& mt) {
                        return mt;
                    }
            };

        template<class T, MajorOrder major_order>
            struct pass_to_expression<const Matrix<T,major_order>&> {
                const Matrix<T,major_order>& 
                    operator()(const Matrix<T,major_order>& mt) {
                        return mt;
                    }
            };

        /*
         * =====================================================================================
         *        Class:  Matrix
         *  Description:  Matrix class
         * =====================================================================================
         */
        template < class T, MajorOrder major_order=RowMajor >
            class Matrix : public MatExpression<Matrix<T,major_order>>
        {
            public:
                /* ====================  LIFECYCLE     ======================================= */

                Matrix (size_t r = 0, size_t c = 0) : _row(r), _col(c), _size(r*c) {    /* constructor */
                    if(_size)
                        _data = new T[_size];
                    else
                        _data = NULL;
#ifndef NDEBUG
                    std::cout << "Constructor Matrix<" << typeid(T).name() 
                        << ", RowMajor>(" << r << ", " << c << " ) at : " << this << std::endl;
#endif
                }                             

                Matrix (T* v, size_t r, size_t c) 
                    : _row(r), _col(c), _size(r*c) {

                        _data = new T[_size];
                        if(is_trivially_copyable<T>::value)
                            std::memcpy(_data,v,_size*sizeof(T));
                        else {
                            for(size_t i =0; i < _size; ++i)
                                _data[i] = v[i];
                        }

#ifndef NDEBUG
                        std::cout << "Constructor Matrix" << typeid(T).name()
                            << ", RowMajor>( std::vector<" << typeid(T).name() <<
                            ">, " << r << ", " << c << " ) at : " << this << std::endl;
#endif
                    }

                template <class E> Matrix(MatExpression<E> const& me) :
                        _row(me.row()), _col(me.col()), _size(me.size()) {

                            if(_size)
                                _data = new T[_size];
                            else
                                _data = NULL;

                            E const& m = me;
                            const auto mj2 = traits<E>::major_order;

                            if(mj2 == RowMajor) {
                                size_t pos;
#ifdef _OPENMP
                                if(_size < MatExpression<E>::ParallelThreshold)
                                    for(pos = 0; pos <  _size; ++pos) {
                                        _data[pos] = m[pos];
                                    }
                                else
#pragma omp parallel for private(pos)
#endif
                                for(pos = 0; pos <  _size; ++pos) {
                                    _data[pos] = m[pos];
                                }

                            } else {
                                size_t pos = (size_t)-1;
                                for(size_t i = 0; i < _col; ++i) {
                                    for(size_t j = 0; j < _row; ++j) {
                                        _data[i+j*_col] = m[++pos];
                                    }
                                }
                            }

#ifndef NDEBUG
                            std::cout << "Constructor Matrix<" << typeid(T).name() <<
                                ", RowMajor> ( " << MatExpression<E>::name() 
                                << " ) at : " << this << std::endl;
#endif
                        }

                Matrix(Matrix<T,major_order> const& m) 
                    : _row(m._row), _col(m._col), _size(m._size) {

                        if(_size)
                            _data = new T[_size];
                        else
                            _data = NULL;

                        if(is_trivially_copyable<T>::value)
                            std::memcpy(_data, m._data, _size*sizeof(T));
                        else {
                            for(size_t i = 0; i < _size; ++i)
                                _data[i] = m._data[i];
                        }
#ifndef NDEBUG
                        std::cout << "Copy Constructor Matrix( " <<
                            Matrix<T,major_order>::name() << " ) at : " 
                            << &m << ", " << this << std::endl;
#endif
                    }
                
                // Move Constructor
                Matrix (Matrix<T,RowMajor>&& m) 
                    : _row(m._row), _col(m._col), _size(m._size), _data(m._data) {
                        m._row = m._col = m._size = 0;
                        m._data = NULL;
#ifndef NDEBUG
                        std::cout << "Move Constructor " << name() << " at " <<
                            &m << ", " << this << std::endl;
#endif
                    }

//                template<class E1, class E2>
//                Matrix (MatExpression<MatMul<E1,E2>&& mm) 


                virtual ~Matrix() {
                    if(_data) {
                        delete[] _data;
                        _data = NULL;
                    }
#ifndef NDEBUG
                    std::cout << "Destructor ~" << name() << "() at " << this << "\n";
#endif
                }


                /*
                 * ===  FUNCTION  ======================================================================
                 *         Name:  Random
                 *  Description: create random matrix
                 * =====================================================================================
                 */

                static Matrix<T,major_order> Random ( size_t r, size_t c, 
                        const T& min, const T& max) {

                    Matrix<T, major_order> mt(r,c);

                    auto ran = boss14420::random::random<T>(min, max);
                    for(size_t i = 0; i < r*c; ++i) {
                        mt[i] = ran.nextRan();
                    }

                    return mt;
                }
                /* -----  end of function Random  ----- */


                /* ====================  ACCESSORS     ======================================= */

                size_t row() const { return _row;}
                size_t col() const { return _col;}
                size_t size() const { return _size;}

                T at(size_t i) const {
                    if(i >= _size) {
                        // throw
                        return 0;
                    }
                    return _data[i];
                }

                T& at(size_t i) {
                    if(i >= _size) {
                        // throw
                    }
                    return _data[i];
                }

                T at(size_t i, size_t j) const {
                    if(i >= _row || j >= _col) {
                        // throw
                        return 0;
                    }
                    return _data[i*_col+j];
                }

                T& at(size_t i, size_t j) {
                    if(i >= _row || j >= _col) {
                        // throw
                    }
                    return _data[i*_col+j];
                }

                T* getRawData() { return _data; }
                const T* getRawData() const { return _data; }


                const Matrix<T,major_order>& eval() const { return *this; }

#ifndef NDEBUG
                static std::string name() {
                    return std::string("Matrix<") + typeid(T).name() 
                        + ", RowMajor>";
                }
#endif


                /* ====================  MUTATORS      ======================================= */

                void resize(size_t r, size_t c) {
                    if(r != 0 && c != 0) {
                        _row = r;
                        _col = c;
                        _size = r*c;
                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                        _data = new T[_size];
                    } else {
                        //throw
                    }
                }

                /* ====================  OPERATORS     ======================================= */

                template<class E>
                    Matrix& operator= (MatExpression<E> const& me) {
                        _row = me.row();
                        _col = me.col();
                        _size = me.size();

                        const auto mj = traits<E>::major_order;

                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                        _data = new T[_size];
                        if(mj == major_order) {
                            size_t i;
#ifdef _OPENMP
                            if(_size < MatExpression<E>::ParallelThreshold)
                                for(i = 0; i < _size; ++i) {
                                    _data[i] = me[i];
                                }
                            else
#pragma omp parallel for private(i)
#endif
                            for(i = 0; i < _size; ++i) {
                                _data[i] = me[i];
                            }
                        } else {
                            size_t pos = (size_t)-1;
                            for(size_t i = 0; i < _col; ++i) {
                                for(size_t j = 0; j < _row; ++j) {
                                    _data[i+j*_col] = me[++pos];
                                }
                            }
                        }

#ifndef NDEBUG
                        std::cout << name() << "::operator=" << std::endl;
#endif
                        return *this;
                    }

                Matrix& operator= (const Matrix<T,RowMajor>& mt) {
                    if(&mt != this) {
                        _row = mt._row;
                        _col = mt._col;
                        _size = mt._size;

                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                        if(_size) {
                            _data = new T[_size];
                        }
                        if(is_trivially_copyable<T>::value)
                            std::memcpy(_data,mt._data,_size*sizeof(T));
                        else {
                            for(size_t i = 0; i < _size; ++i)
                                _data[i] = mt._data[i];
                        }
                    }

                    return *this;
                }
                
                // Move assigment
                Matrix& operator= (Matrix<T,RowMajor>&& mt) {
                    if(&mt != this) {
                        _row = mt._row;
                        _col = mt._col;
                        _size = mt._size;

                        if(_data) {
                            delete[] _data;
                        }

                        _data = mt._data;

                        mt._row = mt._col = mt._size = 0;
                        mt._data = NULL;
                    }
#ifndef NDEBUG
                    std::cout << "Move assignment " << name() << " at " <<
                        &mt << ", " << this << std::endl;
#endif
                    return *this;
                }


                STRONG_INLINE T operator[](size_t i) const { return _data[i];}
                STRONG_INLINE T& operator[](size_t i) { return _data[i];}
                STRONG_INLINE T operator()(size_t i, size_t j) const { // Fortran style m(i,j)
                    return _data[i*_col+j];
                }
                STRONG_INLINE T& operator()(size_t i, size_t j) {
                    return _data[i*_col+j];
                }

                // IO operator
                friend std::istream& operator>> (std::istream& is, Matrix<T,major_order>& mt) {
                    if(&is == &std::cin || !is_trivially_copyable<T>::value) { // Standard input
                        std::cout << "Type matrix's size (row & col): ";
                        is >> mt._row >> mt._col;
                        mt._size = mt._row * mt._col;
                        if(mt._data) {
                            delete[] mt._data;
                            mt._data = NULL;
                        }
                        mt._data = new T[mt._size];

                        std::istream_iterator<T> ii_begin(is), ii_end;
                        //                T*::iterator vi = mt._data.begin();
                        auto vi = mt._data;
                        for(size_t i = 0; i < mt._row; ++i) {
                            std::cout << "Matrix's row " << i << " : ";
                            std::copy(ii_begin, ii_end, vi);
                            vi += mt._col;
                        }
                    } else {
                        char mj;
                        is.read((char*)&mt._row, sizeof(size_t))
                            .read((char*)&mt._col, sizeof(size_t))
                            .read(&mj,1);
                        mt._size = mt._row * mt._col;
                        if(mt._data) {
                            delete[] mt._data;
                            mt._data = NULL;
                        }
                        mt._data = new T[mt._size];
                        if(mj == RowMajor) {
                            for(size_t i = 0; i < mt._size; ++i) {
                                is.read((char*)&mt._data[i], sizeof(T));
                            }
                        } else {
                            for(size_t i = 0; i < mt._col; ++i) {
                                for(size_t j = 0; j < mt._row; ++j) {
                                    is.read((char*)&mt._data[i+j*mt._col],sizeof(T));
                                }
                            }
                        }
                    }

                    return is;
                }

            protected:
                /* ====================  DATA MEMBERS  ======================================= */

            private:
                /* ====================  DATA MEMBERS  ======================================= */

                size_t _row, _col, _size;
                T* _data;

        }; /* ----------  end of template class Matrix  ---------- */



        /*
         * =====================================================================================
         *        Class:  Matrix
         *  Description:  template specialization for ColumnMajor
         *  Matrix
         * =====================================================================================
         */
        template < class T >
            class Matrix<T, ColumnMajor> : public MatExpression<Matrix<T,ColumnMajor>> {
            public:
                /* ====================  LIFECYCLE     ======================================= */

                Matrix (size_t r = 0, size_t c = 0) : _row(r), _col(c), _size(r*c) {    /* constructor */
                    if(_size)
                        _data = new T[_size];
                    else
                        _data = NULL;
#ifndef NDEBUG
                    std::cout << "Constructor Matrix<" << typeid(T).name() 
                        << ", ColumnMajor>(" << r << ", " << c << " ) at : "
                        << this << std::endl;
#endif
                }                             

                Matrix (T* v, size_t r, size_t c) 
                    : _row(r), _col(c), _size(r*c) {

                        _data = new T[_size];
                        if(is_trivially_copyable<T>::value)
                            std::memcpy(_data, v, _size*sizeof(T));
                        else {
                            for(size_t i = 0; i < _size; ++i)
                                _data[i] = v[i];
                        }
#ifndef NDEBUG
                        std::cout << "Constructor Matrix ( std::vector<" << typeid(T).name() <<
                            "ColumnMajor>, " << r << ", " << c << " ) at : " 
                            << this << std::endl;
#endif
                    }

                template <class E>
                    Matrix(MatExpression<E> const& me) : 
                        _row(me.row()), _col(me.col()), _size(me.size()) {

                            if(_size)
                                _data = new T[_size];
                            else
                                _data = NULL;

                            E const& m = me;
                            auto const mj2 = traits<E>::major_order;

                            if(mj2 == ColumnMajor) {
                                size_t pos;
#ifdef _OPENMP
                                if(_size < MatExpression<E>::ParallelThreshold)
                                    for(pos = 0; pos < _size; ++pos) {
                                        _data[pos] = m[pos];
                                    }
                                else
#pragma omp parallel for private(pos)
#endif
                                for(pos = 0; pos < _size; ++pos) {
                                    _data[pos] = m[pos];
                                }

                            } else {
                                size_t pos = (size_t)-1;
                                for(size_t i = 0; i < _row; ++i) {
                                    for(size_t j = 0; j < _col; ++j) {
                                        _data[i+j*_row] = m[++pos];
                                    }
                                }
                            }

#ifndef NDEBUG
                            std::cout << "Constructor Matrix<" << typeid(T).name()
                                << ", ColumnMajor> ( " << MatExpression<E>::name()
                                << " ) at : " << this << std::endl;
#endif
                        }

                Matrix(Matrix<T,ColumnMajor> const& m) 
                    : _row(m._row), _col(m._col), _size(m._size) {

                        if(_size)
                            _data = new T[_size];
                        else
                            _data = NULL;
                        if(is_trivially_copyable<T>::value)
                            std::memcpy(_data, m._data, _size*sizeof(T));
                        else {
                            for(size_t i = 0; i < _size; ++i)
                                _data[i] = m._data[i];
                        }
#ifndef NDEBUG
                        std::cout << "Copy Constructor Matrix<" << typeid(T).name()
                            << ", ColumnMajor> ( " << Matrix<T,ColumnMajor>::name() 
                            << " ) at : " << &m << ", " << this << std::endl;
#endif
                    }

                // Move Constructor
                Matrix (Matrix<T,ColumnMajor>&& m) 
                    : _row(m._row), _col(m._col), _size(m._size), _data(m._data) {
                        m._row = m._col = m._size = 0;
                        m._data = NULL;
#ifndef NDEBUG
                        std::cout << "Move Constructor " << name() << " at " <<
                            &m << ", " << this << std::endl;
#endif
                    }

                virtual ~Matrix() {
                    if(_data) {
                        delete[] _data;
                        _data = NULL;
                    }
#ifndef NDEBUG
                    std::cout << "Destructor ~" << name() << "() at " << this << "\n";
#endif
                }


                /*
                 * ===  FUNCTION  ======================================================================
                 *         Name:  Random
                 *  Description: create random matrix
                 * =====================================================================================
                 */

                static Matrix<T,ColumnMajor> Random ( size_t r, size_t c, 
                        const T& min, const T& max) {

                    Matrix<T, ColumnMajor> mt(r,c);

                    auto ran = boss14420::random::random<T>(min, max);
                    for(size_t i = 0; i < r*c; ++i) {
                        mt[i] = ran.nextRan();
                    }

                    return mt;
                }
                /* -----  end of function Random  ----- */


                /* ====================  ACCESSORS     ======================================= */

                size_t row() const { return _row;}
                size_t col() const { return _col;}
                size_t size() const { return _size;}

                T at(size_t i) const {
                    if(i >= _size) {
                        // throw
                        return 0;
                    }
                    return _data[i];
                }

                T& at(size_t i) {
                    if(i >= _size) {
                        // throw
                    }
                    return _data[i];
                }

                T at(size_t i, size_t j) const {
                    if(i >= _row || j >= _col) {
                        // throw
                        return 0;
                    }
                    return _data[i+_row*j];
                }

                T& at(size_t i, size_t j) {
                    if(i >= _row || j >= _col) {
                        // throw
                    }
                    return _data[i+_row*j];
                }

                T* getRawData() { return _data; }
                const T* getRawData() const { return _data; }

                const Matrix<T,ColumnMajor>& eval() const { return *this; }

#ifndef NDEBUG
                static std::string name() {
                    return std::string("Matrix<") + typeid(T).name() 
                        + ", ColumnMajor>";
                }
#endif


                /* ====================  MUTATORS      ======================================= */

                void resize(size_t r, size_t c) {
                    if(r != 0 && c != 0) {
                        _row = r;
                        _col = c;
                        _size = r*c;
                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                        _data = new T[_size];
                    } else {
                        //throw
                    }
                }

                /* ====================  OPERATORS     ======================================= */

                template<class E>
                    Matrix& operator= (MatExpression<E> const& me) {
                        _row = me.row();
                        _col = me.col();
                        _size = me.size();

                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                        _data = new T[_size];

                        auto const mj = traits<E>::major_order;

                        if(mj == ColumnMajor) {
                            size_t i;
#ifdef _OPENMP
                            if(_size < MatExpression<E>::ParallelThreshold)
                                for(i = 0; i < _size; ++i) {
                                    _data[i] = me[i];
                                }
                            else
#pragma omp parallel for private(i)
#endif
                            for(i = 0; i < _size; ++i) {
                                _data[i] = me[i];
                            }
                        } else {
                            size_t pos = -1;
                            for(size_t i = 0; i < _row; ++i) {
                                for(size_t j = 0; j < _col; ++j) {
                                    _data[i+j*_row] = me[++pos];
                                }
                            }
                        }
#ifndef NDEBUG
                        std::cout << name() << "::operator=" << std::endl;
#endif

                        return *this;
                    }

                Matrix& operator= (const Matrix<T,ColumnMajor>& mt) {
                    if(&mt != this) {
                        _row = mt._row;
                        _col = mt._col;
                        _size = mt._size;

                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                        if(_size) {
                            _data = new T[_size];
                        }
                        if(is_trivially_copyable<T>::value)
                            std::memcpy(_data,mt._data,_size*sizeof(T));
                        else {
                            for(size_t i = 0; i < _size; ++i)
                                _data[i] = mt._data[i];
                        }
                    }

                    return *this;
                }

                // Move assigment
                Matrix& operator= (Matrix<T,ColumnMajor>&& mt) {
                    if(&mt != this) {
                        std::swap(*this,mt);
                    }
#ifndef NDEBUG
                    std::cout << "Move assignment " << name() << " at " <<
                        &mt << ", " << this << std::endl;
#endif
                    return *this;
                }

                STRONG_INLINE T operator[](size_t i) const { return _data[i];}
                STRONG_INLINE T& operator[](size_t i) { return _data[i];}
                STRONG_INLINE T operator()(size_t i, size_t j) const { // Fortran style m(i,j)
                    return _data[i+j*_row];
                }
                STRONG_INLINE T& operator()(size_t i, size_t j) {
                    return _data[i+j*_row];
                }

                // IO operator

                
                friend std::istream& operator>> (std::istream& is, Matrix<T,ColumnMajor>& mt) {
                    if(&is == &std::cin || !is_trivially_copyable<T>::value) { // Standard input
                        std::cout << "Type matrix's size (row & col): ";
                        is >> mt._row >> mt._col;
                        mt._size = mt._row * mt._col;
                        mt.resize(mt._row,mt._col);

                        for(size_t i = 0; i < mt._row; ++i) {
                            for(size_t j = 0; j < mt._col; ++j) {
                                is >> mt(i,j);
                            }
                        }
                    } else {
                        char mj;
                        is.read((char*)&mt._row, sizeof(size_t))
                            .read((char*)&mt._col, sizeof(size_t))
                            .read(&mj,1);

                        mt._size = mt._row * mt._col;
                        if(mt._data) {
                            delete[] mt._data;
                            mt._data = NULL;
                        }
                        mt._data = new T[mt._size];

                        if(mj == ColumnMajor) {
                            for(size_t i = 0; i < mt._size; ++i) {
                                is.read((char*)&mt._data[i], sizeof(T));
                            }
                        } else {
                            for(size_t i = 0; i < mt._row; ++i) {
                                for(size_t j = 0; j < mt._col; ++j) 
                                    is.read((char*)&mt._data[i+j*mt._row], sizeof(T));
                            }
                        }
                    }

                    return is;
                }

            protected:
                /* ====================  DATA MEMBERS  ======================================= */

            private:
                /* ====================  DATA MEMBERS  ======================================= */

                size_t _row, _col, _size;
                T* _data;


        }; /* ----------  end of template class Matrix  ---------- */



    }
}
#endif
