/*
 * =====================================================================================
 *
 *       Filename:  vector.hpp
 *
 *    Description:  Vector
 *
 *        Version:  1.0
 *        Created:  11/27/2011 03:44:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  BOSS14420 (boss14420), boss14420@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */
#include <vector>
#include <string>

#ifdef DEBUG
    #include <iostream>
    #include <typeinfo>
#endif

#include <cassert>
#include "macro.h"
 
template <typename T, typename E>
// A CRTP base class for Vecs with a size and indexing:
class VecExpression {
    public:
        typedef typename std::vector<T>         container_type;
        typedef typename container_type::size_type   size_type;
        typedef typename container_type::value_type  value_type;
        typedef typename container_type::reference   reference;

        size_type size() const { return ((E const&)*this).size(); }
        inline value_type operator[](size_type i) const { return ((E const&)*this)[i]; }

#ifdef DEBUG
        virtual std::string name() const {
            return std::string("VecExpression< ") + typeid(T).name() 
                + ", " + ((E const&)*this).name() + " >";
        }
#endif

        operator E&() { return reinterpret_cast<E&>(*this); }
        operator E const&() const { return reinterpret_cast<E const&>(*this); }
};

// The actual Vec class:
template<typename T>
class Vec : public VecExpression<T, Vec<T> > {
    typedef typename VecExpression<T, Vec<T> >::container_type container_type;
    container_type _data;

public:
    typedef typename VecExpression<T, Vec<T> >::size_type size_type;
    typedef typename VecExpression<T, Vec<T> >::value_type value_type;
    typedef typename VecExpression<T, Vec<T> >::reference reference;

#ifdef DEBUG
    virtual std::string name() const {
        return std::string("Vec< ") + typeid(T).name() + " >";
    }
#endif

    inline reference operator[](size_type i) { return _data[i]; }
    inline value_type operator[](size_type i) const { return _data[i]; }
    size_type size() const { return _data.size(); }

    Vec(size_type n) : _data(n) {} // Construct a given size:

    // Construct from any VecExpression:
    template <typename E>
        Vec(VecExpression<T, E> const& vec) {
            E const& v = vec;
            _data.resize(v.size());
            for (size_type i = 0; i != v.size(); ++i) {
                _data[i] = v[i];
            }
#ifdef DEBUG
            std::cout << "Call constructor Vec(VecExpression< " <<
                typeid(T).name() << ", " << v.name() << " >)\n";
#endif
        }

    template <typename E>
        Vec& operator= (VecExpression<T, E> const& vee) {
            _data.resize(vee.size());
            for(size_type i = 0; i < vee.size(); ++i) {
                _data[i] = vee[i];
            }
        }

    Vec(container_type &dt) {
        _data = dt;
#ifdef DEBUG
        std::cout << "Call constructor Vec(" << 
            typeid(container_type).name() << ")\n";
#endif
    }
};

template <typename T, typename E1, typename E2>
class VecDifference : public VecExpression<T, VecDifference<T, E1, E2> > {
    typedef typename Vec<T>::size_type size_type;
    typedef typename Vec<T>::value_type value_type;

    E1 const& _u;
    E2 const& _v;
    
public:
    VecDifference(VecExpression<T, E1> const& u, VecExpression<T, E2> const& v) : _u(u), _v(v) {
        assert(u.size() == v.size());

#ifdef DEBUG
        std::cout << "Call constructor VecDifference(VecExpression< " << typeid(T).name() <<
            ", " << _u.name() << " >, VecExpression< " << typeid(T).name() << ", " <<
            _v.name() << " >)\n";
#endif
    }
    size_type size() const { return _v.size(); }
    inline value_type operator[](size_type i) const { return _u[i] - _v[i]; }

#ifdef DEBUG
    virtual std::string name()  const {
        return std::string("VecDifference< ") + typeid(T).name() +
            ", " + _u.name() + ", " + _v.name() + " >";
    }
#endif
};

template <typename T, typename E1, typename E2>
class VecAdd : public VecExpression<T, VecAdd<T, E1, E2> > {
    typedef typename Vec<T>::size_type size_type;
    typedef typename Vec<T>::value_type value_type;

    E1 const& _u;
    E2 const& _v;
    
public:
    VecAdd(VecExpression<T, E1> const& u, VecExpression<T, E2> const& v) : _u(u), _v(v) {
        assert(u.size() == v.size());

#ifdef DEBUG
        std::cout << "Call constructor VecAdd(VecExpression< " << typeid(T).name() <<
            ", " << _u.name() << " >, VecExpression< " << typeid(T).name() << ", " <<
            _v.name() << " >)\n";
#endif
    }
    size_type size() const { return _v.size(); }
    inline value_type operator[](size_type i) const { return _u[i] + _v[i]; }

#ifdef DEBUG
    virtual std::string name() const {
        return std::string("VecAdd< ") + typeid(T).name() +
            ", " + _u.name() + ", " + _v.name() + " >";
    }
#endif
};

template <typename T, typename E>
class VecScaled : public VecExpression<T, VecScaled<T, E> > {
    typedef typename Vec<T>::size_type size_type;
    typedef typename Vec<T>::value_type value_type;

    T _alpha; 
    E const& _v;

public:
    VecScaled(T alpha, VecExpression<T,E> const& v) : _alpha(alpha), _v(v) {
#ifdef DEBUG
        std::cout << "Call constructor VecScaled(" << typeid(T).name() << 
            ", VecExpression< " << _v.name() << " >)\n";
#endif
    }
    size_type size() const { return _v.size(); }
    inline value_type operator[](size_type i) const { return _alpha * _v[i]; }

#ifdef DEBUG
    virtual std::string name() const {
        return std::string("VecScaled< ") + typeid(T).name() + ", " + _v.name() + " >";
    }
#endif
};

// Now we can overload operators:

template <typename T, typename E1, typename E2>
VecDifference<T,E1,E2> const
operator-(VecExpression<T, E1> const& u, VecExpression<T, E2> const& v) {
    return VecDifference<T,E1,E2>(u,v);
}

template <typename T, typename E1, typename E2>
VecAdd<T,E1,E2> const
operator+(VecExpression<T, E1> const& u, VecExpression<T, E2> const& v) {
    return VecAdd<T,E1,E2>(u,v);
}

template <typename T, typename E>
VecScaled<T,E> const
operator*(T alpha, VecExpression<T,E> const& v) {
    return VecScaled<T,E>(alpha,v);
}

template <typename T, typename E>
VecScaled<T,E> const
operator*(VecExpression<T,E> const& v, T alpha) {
    return VecScaled<T,E>(alpha,v);
}

template <typename T, typename E1, typename E2>
T operator*(VecExpression<T,E1> const& u, VecExpression<T, E2> const& v) {
    assert(u.size() == v.size());
    typename Vec<T>::size_type sz = u.size(), i = 0;
    T res = 0;
    for(; i < sz; ++i) {
        res += u[i] * v[i];
    }

    return res;
}
