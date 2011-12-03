#ifndef _TRIVIAL_MATRIX_

#include <cstring>
#include <cstdlib>

#include "../classes/random.hpp"

namespace boss14420 {
    namespace trivial {
        template<class T>
            class Matrix {
                public:
                    Matrix (const Matrix<T>& mt) 
                        : _row(mt._row), _col(mt._col), _size(mt._size) {
                            if(_size)
                                _data = new T[_size];
                            else
                                _data = NULL;

                            std::memcpy(_data, mt._data, _size*sizeof(T));
                        }

                    Matrix (size_t r = 0, size_t c = 0) : _row(r), _col(c), _size(r*c) {
                        if(_size)
                            _data = new T[_size];
                        else
                            _data = NULL;
                    }

                    virtual ~Matrix() {
                        if(_data) {
                            delete[] _data;
                            _data = NULL;
                        }
                    }

                    Matrix<T>& operator=(const Matrix<T>& mt) {
                        _row = mt._row;
                        _col = mt._col;
                        _size = mt._size;

                        if(_data)
                            delete[] _data;
                        if(_size)
                            _data = new T[_size];
                        else
                            _data = NULL;

                        std::memcpy(_data, mt._data, _size*sizeof(T));

                        return *this;
                    }

                    T operator[](size_t i) const { return _data[i];}
                    T& operator[](size_t i) { return _data[i];}

                    T operator()(size_t i, size_t j) const {
                        return _data[i*_col+j];
                    }
                    T& operator()(size_t i, size_t j) {
                        return _data[i*_col+j];
                    }

                    Matrix<T> operator+(const Matrix<T>& mt) const {
                        Matrix<T> mres(_row, _col);

                        if(mt._row == _row && mt._col == _col) {
                            T* _data2 = mt._data;
                            T* _datar = mres._data;
                            for(size_t i = 0; i < _size; ++i)
                                _datar[i] = _data[i] + _data2[i];
                        }

                        return mres;
                    }

                    static Matrix<T> Random ( size_t r, size_t c, 
                        const T& min, const T& max) {

                        Matrix<T> mt(r,c);

                        auto ran = boss14420::random::random<T>(min, max);
                        for(size_t i = 0; i < r*c; ++i) {
                            mt[i] = ran.nextRan();
                        }

                        return mt;
                    }

                private:
                    T* _data;
                    size_t _row, _col, _size;
            };
    }
}

#endif
