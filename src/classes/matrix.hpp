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
#include <cstdlib>
#include <ctime>
#include <cstring>

#include "random.hpp"


namespace boss14420 {
    namespace matrix {

        enum MajorOrder { RowMajor = 0, ColumnMajor = 1 };

        template<class T, MajorOrder mj> 
            class Matrix;

        /*
         * =====================================================================================
         *        Class:  MatExpression
         *  Description:  template Matrix expression class
         * =====================================================================================
         */
        template < class T, MajorOrder major_order, class E >
            class MatExpression
            {
                public:
                    /* ====================  LIFECYCLE     ======================================= */

                    /* ====================  ACCESSORS     ======================================= */
                    size_t row() const {
                        return static_cast<E const&>(*this).row();
                    }

                    size_t col() const {
                        return static_cast<E const&>(*this).col();
                    }

                    size_t size() const {
                        return static_cast<E const&>(*this).size();
                    }

                    T at(size_t i) const {
                        return static_cast<E const&>(*this).at(i);
                    }

                    T at(size_t i, size_t j) const {
                        return static_cast<E const&>(*this).at(i,j);
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

                    T operator[](size_t i) const {
                        return static_cast<E const&>(*this)[i];
                    }

                    T operator()(size_t i, size_t j) const {
                        return static_cast<E const&>(*this)(i, j);
                    }

                    operator E&() {
                        return static_cast<E&>(*this);
                    }

                    operator E const&() const {
                        return static_cast<E const&>(*this);
                    }

                    friend std::ostream& operator<< (std::ostream& os, const MatExpression<T,major_order,E>& me) {
                        //            if(typeid(E) != typeid(Matrix<T>)) {
                        //                return os << static_cast<E&>(me).eval();
                        //            } else {
                        //                return os << static_cast<Matrix<T> const&>(me);
                        //            }
                        return os << static_cast<Matrix<T,major_order> const&>(me);
                    }

                protected:
                    /* ====================  DATA MEMBERS  ======================================= */

                private:
                    /* ====================  DATA MEMBERS  ======================================= */

            }; /* ----------  end of template class MatExpression  ---------- */


        /*
         * =====================================================================================
         *        Class:  Matrix
         *  Description:  Matrix class
         * =====================================================================================
         */
        template < class T, MajorOrder major_order=RowMajor >
            class Matrix : public MatExpression<T,major_order,Matrix<T> >
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
                        std::memcpy(_data,v,_size*sizeof(T));

#ifndef NDEBUG
                        std::cout << "Constructor Matrix" << typeid(T).name()
                            << ", RowMajor>( std::vector<" << typeid(T).name() <<
                            ">, " << r << ", " << c << " ) at : " << this << std::endl;
#endif
                    }

                template <class T2, MajorOrder mj2, class E>
                    Matrix(MatExpression<T2, mj2, E> const& me) :
                        _row(me.row()), _col(me.col()), _size(me.size()) {

                            if(_size)
                                _data = new T[_size];
                            else
                                _data = NULL;

                            E const& m = me;
                            if(mj2 == RowMajor) {
                                for(size_t pos = 0; pos != _size; ++pos) {
                                    _data[pos] = m[pos];
                                }

                            } else {
                                size_t pos = -1;
                                for(size_t i = 0; i < _col; ++i) {
                                    for(size_t j = 0; j < _row; ++j) {
                                        _data[i+j*_col] = m[++pos];
                                    }
                                }
                            }

#ifndef NDEBUG
                            std::cout << "Constructor Matrix<" << typeid(T).name() <<
                                ", RowMajor> ( " << MatExpression<T2,major_order,E>::name() 
                                << " ) at : " << this << std::endl;
#endif
                        }

                Matrix(Matrix<T,major_order> const& m) 
                    : _row(m._row), _col(m._col), _size(m._size) {

                        if(_size)
                            _data = new T[_size];
                        else
                            _data = NULL;
                        std::memcpy(_data, m._data, _size*sizeof(T));
#ifndef NDEBUG
                        std::cout << "Copy Constructor Matrix( " <<
                            Matrix<T,major_order>::name() << " ) at : " 
                            << &m << ", " << this << std::endl;
#endif
                    }

                virtual ~Matrix() {
                    if(_data) {
                        delete[] _data;
                        _data = NULL;
                    }
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

                template<class T2, MajorOrder mj, class E>
                    Matrix& operator= (MatExpression<T2, mj, E> const& me) {
                        _row = me.row();
                        _col = me.col();
                        _size = me.size();

                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                        _data = new T[_size];
                        if(mj == major_order) {
                            for(size_t i = 0; i < _size; ++i) {
                                _data[i] = me[i];
                            }
                        } else {
                            size_t pos = -1;
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
                    std::memcpy(_data,mt._data,_size*sizeof(T));

                    return *this;
                }

                T operator[](size_t i) const { return _data[i];}
                T& operator[](size_t i) { return _data[i];}
                T operator()(size_t i, size_t j) const { // Fortran style m(i,j)
                    return _data[i*_col+j];
                }
                T& operator()(size_t i, size_t j) {
                    return _data[i*_col+j];
                }

                // IO operator

                friend std::ostream& operator<< (std::ostream& os, const Matrix<T,major_order>& mt) {
                    if(&os == &std::cout) { // standard input
                        os << mt._row << "x" << mt._col << std::endl;

                        //                typename T*::const_iterator ii = mt._data.begin();
                        auto ii = mt._data;
                        for(size_t i = 0; i < mt._row; ++i) {
                            std::copy(ii, ii+mt._col, std::ostream_iterator<T>(os, " "));
                            ii+=mt._col;
                            os << std::endl;
                        }
                    } else {
                        char mj = major_order;
                        os.write((char*)&mt._row, sizeof(size_t))
                            .write((char*)&mt._col, sizeof(size_t))
                            .write(&mj, 1);

                        for(size_t i = 0; i < mt._size; ++i) {
                            os.write((char*)&mt._data[i], sizeof(T));
                        }
                    }

                    return os;
                }

                friend std::istream& operator>> (std::istream& is, Matrix<T,major_order>& mt) {
                    if(&is == &std::cin) { // Standard input
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
            class Matrix<T, ColumnMajor> : MatExpression<T, ColumnMajor, Matrix<T,ColumnMajor>> {
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
                        std::memcpy(_data, v, _size*sizeof(T));                            
#ifndef NDEBUG
                        std::cout << "Constructor Matrix ( std::vector<" << typeid(T).name() <<
                            "ColumnMajor>, " << r << ", " << c << " ) at : " 
                            << this << std::endl;
#endif
                    }

                template <class T2, MajorOrder mj2, class E>
                    Matrix(MatExpression<T2, mj2, E> const& me) : 
                        _row(me.row()), _col(me.col()), _size(me.size()) {

                            if(_size)
                                _data = new T[_size];
                            else
                                _data = NULL;

                            E const& m = me;
                            if(mj2 == ColumnMajor) {
                                for(size_t pos = 0; pos != _size; ++pos) {
                                    _data[pos] = m[pos];
                                }

                            } else {
                                size_t pos = -1;
                                for(size_t i = 0; i < _row; ++i) {
                                    for(size_t j = 0; j < _col; ++j) {
                                        _data[i+j*_row] = m[++pos];
                                    }
                                }
                            }

#ifndef NDEBUG
                            std::cout << "Constructor Matrix<" << typeid(T).name()
                                << ", ColumnMajor> ( " << MatExpression<T2,ColumnMajor,E>::name()
                                << " ) at : " << this << std::endl;
#endif
                        }

                Matrix(Matrix<T,ColumnMajor> const& m) 
                    : _row(m._row), _col(m._col), _size(m._size) {

                        if(_size)
                            _data = new T[_size];
                        else
                            _data = NULL;
                        std::memcpy(_data, m._data, _size*sizeof(T));
#ifndef NDEBUG
                        std::cout << "Copy Constructor Matrix<" << typeid(T).name()
                            << ", ColumnMajor> ( " << Matrix<T,ColumnMajor>::name() 
                            << " ) at : " << &m << ", " << this << std::endl;
#endif
                    }

                virtual ~Matrix() {
                    if(_data) {
                        delete[] _data;
                        _data = NULL;
                    }
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

                template<class T2, MajorOrder mj, class E>
                    Matrix& operator= (MatExpression<T2, mj, E> const& me) {
                        _row = me.row();
                        _col = me.col();
                        _size = me.size();

                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                        _data = new T[_size];

                        if(mj == ColumnMajor) {
                            for(size_t i = 0; i < _size; ++i) {
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
                    std::memcpy(_data,mt._data,_size*sizeof(T));

                    return *this;
                }

                T operator[](size_t i) const { return _data[i];}
                T& operator[](size_t i) { return _data[i];}
                T operator()(size_t i, size_t j) const { // Fortran style m(i,j)
                    return _data[i+j*_row];
                }
                T& operator()(size_t i, size_t j) {
                    return _data[i+j*_row];
                }

                // IO operator

                friend std::ostream& operator<< (std::ostream& os, const Matrix<T,ColumnMajor>& mt) {
                    if(&os == &std::cout) { // standard input
                        os << mt._row << "x" << mt._col << std::endl;

                        //                typename T*::const_iterator ii = mt._data.begin();
                        auto ii = mt._data;
                        for(size_t i = 0; i < mt._row; ++i) {
                            for(size_t j = 0; j < mt._col; ++j) {
                                os << mt(i,j) << " ";
                            }
                            os << std::endl;
                        }
                    } else {
                        char mj = (char)ColumnMajor;
                        os.write((char*)&mt._row, sizeof(size_t))
                            .write((char*)&mt._col, sizeof(size_t))
                            .write((char*)&mj, sizeof(MajorOrder))
                            .write(&mj,1);

                        for(size_t i = 0; i < mt._size; ++i) {
                            os.write((char*)&mt._data[i], sizeof(T));
                        }
                        //                std::ostreambuf_iterator<T> obi(os);
                        //                std::copy(mt._data.begin(), mt._data.end(), obi);
                    }

                    return os;
                }

                friend std::istream& operator>> (std::istream& is, Matrix<T,ColumnMajor>& mt) {
                    if(&is == &std::cin) { // Standard input
                        std::cout << "Type matrix's size (row & col): ";
                        is >> mt._row >> mt._col;
                        mt._size = mt._row * mt._col;
                        mt._data.resize(mt._size);
                        if(mt._data) {
                            delete[] mt._data;
                            mt._data = NULL;
                        }
                        mt._data = new T[mt._size];

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