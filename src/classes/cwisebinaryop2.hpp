/*
 * =====================================================================================
 *
 *       Filename:  CWiseBinaryOp2.hpp
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
 *RowMajor>(3,
 * =====================================================================================
 */

#ifndef _CWISEBINARYOP2_HPP_
#define _CWISEBINARYOP2_HPP_

#include "matrix.hpp"
#include "functor.hpp"

namespace boss14420 {
namespace matrix {

//template<class T1, class T2>
//    auto _op(T1 t1, T2 t2) -> decltype(t1+t2) {
//        return t1+t2;
//    }

template<class _E1, class _E2>
struct traits<CWiseBinaryOp2< _E1, _E2>> {
//	typedef _BinaryOp BinaryOp;
	typedef _E1 E1;
	typedef _E2 E2;

	typedef typename traits<E1>::T T1;
	typedef typename traits<E2>::T T2;

        typedef scalar_sum_op<T1,T2> BinaryOp;
        typedef typename BinaryOp::T T;
//	typedef decltype(BinaryOp()(T1(),T2())) T;
//        typedef decltype(_op(T1(),T2())) T;

	static const MajorOrder major_order = traits<E1>::major_order;
};

/*
 * =====================================================================================
 *        Class:  CWiseBinaryOp2
 *  Description:  Matrix Add Expression class
 * =====================================================================================
 */
template<class E1, class E2>
class CWiseBinaryOp2: public MatExpression<
		CWiseBinaryOp2<E1, E2>> {
public:
	typedef typename traits<E1>::T T1;
	typedef typename traits<E2>::T T2;
        typedef scalar_sum_op<T1,T2> BinaryOp;

	static const MajorOrder major_order = traits<E1>::major_order;
	typedef decltype(BinaryOp()(T1(),T2())) T;
//        typedef decltype(_op(T1(),T2())) T;

	/* ====================  LIFECYCLE     ======================================= */
	CWiseBinaryOp2(MatExpression<E1> const& u, MatExpression<E2> const& v) :
			_u(u), _v(v), _tmpres(), _evaluted(false)/* , _op(BinaryOp())*/ { /* constructor */
#ifndef NDEBUG
		std::cout << "Constructor CWiseBinaryOp2 ( " << MatExpression<E1>::name()
				<< ", " << MatExpression<E2>::name() << " ) at " << this
				<< " \n";
		std::cout << "&u = " << &u << ", &v = " << &v << "\n";
		std::cout << std::endl;
//                          std::cout << "u = " << u << "\nv = " << v << "\n\n";
#endif

	}

	CWiseBinaryOp2(const CWiseBinaryOp2<E1, E2>& mp) :
			_u(mp._u), _v(mp._v), _tmpres(mp._tmpres), 
                        _evaluted(mp._evaluted)/* , _op( mp._op) */{
#ifndef NDEBUG
		std::cout << "Copy Constructor CWiseBinaryOp2 ( " << name() << " )\n\n";
#endif
	}

	/* ====================  ACCESSORS     ======================================= */

	size_t row() const {
		return _u.row();
	}
	size_t col() const {
		return _u.col();
	}
	size_t size() const {
		return _u.size();
	}

	STRONG_INLINE T at(size_t i) const {
		return _op(_u.at(i), _v.at(i));
	}

	STRONG_INLINE T at(size_t i, size_t j) const {
		return _op(_u.at(i), _v.at(i));
	}

	const Matrix<T, major_order>& eval() {
		if (!_evaluted) {
			_tmpres.resize(row(), col());
                        size_t i;
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
			for (i = 0; i < size(); ++i)
				_tmpres[i] = _op(_u[i], _v[i]);

		}
#ifndef NDEBUG
		std::cout << name() << " at " << this << std::endl;
		std::cout << "_tmpres = " << _tmpres << std::endl << std::endl;
#endif
		return _tmpres;
	}

#ifndef NDEBUG
	static std::string name() {
		return std::string("CWiseBinaryOp2<") + BinaryOp::name() + ", "
				+ E1::name() + ", " + E2::name() + ">";
	}
#endif

	/* ====================  MUTATORS      ======================================= */

	/* ====================  OPERATORS     ======================================= */

	STRONG_INLINE T operator[](size_t i) const {
		return _op(_u[i], _v[i]);
	}

	STRONG_INLINE T operator()(size_t i, size_t j) const {
		return _op(_u(i, j), _v(i, j));
	}

protected:
	/* ====================  DATA MEMBERS  ======================================= */

private:
	/* ====================  DATA MEMBERS  ======================================= */
	E1 const& _u;
	E2 const& _v;
	Matrix<T, major_order> _tmpres;
	bool _evaluted;

	static const BinaryOp _op;

};
/* ----------  end of template class CWiseBinaryOp2  ---------- */

/*
 * =====================================================================================
 *        Class:  CWiseBinaryOp2
 *  Description:  Matrix Add Expression class
 * =====================================================================================
 */
//template<class BinaryOp, class E1, class E2>
//class CWiseBinaryOp2<BinaryOp, E1, E2, true> : public MatExpression<
//		CWiseBinaryOp2<BinaryOp, E1, E2, true>> {
//public:
//	typedef typename traits<E1>::T T1;
//	typedef typename traits<E2>::T T2;
//
//	static const MajorOrder major_order = traits<E1>::major_order;
//	typedef decltype(BinaryOp()(T1(),T2())) T;
//
//	/* ====================  LIFECYCLE     ======================================= */
//
//	CWiseBinaryOp2(MatExpression<E1> const& u, MatExpression<E2> const& v) :
//			_u(u), _v(v), _tmpres(), _evaluted(false)/* , _op(BinaryOp())*/ {
//#ifndef NDEBUG
//		std::cout << "Constructor CWiseBinaryOp2 ( " << MatExpression<E1>::name()
//				<< ", " << MatExpression<E2>::name() << ", bool ) at " << this
//				<< " \n";
//		std::cout << "&u = " << &u << ", &v = " << &v << "\n";
//		std::cout << "&_u = " << &_u << ", &_v = " << &_v << "\n";
//		std::cout << std::endl;
//		//                          std::cout << "u = " << u << "\nv = " << v << "\n\n";
//#endif
//	}
//
//	CWiseBinaryOp2(const CWiseBinaryOp2<BinaryOp, E1, E2, true>& mp) :
//			_u(mp._u), _v(mp._v), _tmpres(mp._tmpres), 
//                        _evaluted(mp._evaluted)/*, _op( mp._op)*/ {
//#ifndef NDEBUG
//		std::cout << "Copy Constructor CWiseBinaryOp2 ( " << name() << " )\n\n";
//#endif
//	}
//
//	/* ====================  ACCESSORS     ======================================= */
//
//	size_t row() const {
//		return _u.row();
//	}
//	size_t col() const {
//		return _u.col();
//	}
//	size_t size() const {
//		return _u.size();
//	}
//
//	STRONG_INLINE T at(size_t i) const {
//		return _op(_u.at(i), _v.at(i));
//	}
//
//	STRONG_INLINE T at(size_t i, size_t j) const {
//		return _op(_u.at(i), _v.at(i));
//	}
//
//	const Matrix<T, major_order>& eval() {
//		if (!_evaluted) {
//			_tmpres.resize(row(), col());
//			for (size_t i = 0; i < size(); ++i)
//				_tmpres[i] = _op(_u[i], _v[i]);
//
//		}
//#ifndef NDEBUG
//		std::cout << name() << " at " << this << std::endl;
//		std::cout << "_tmpres = " << _tmpres << std::endl << std::endl;
//#endif
//		return _tmpres;
//	}
//
//#ifndef NDEBUG
//	static std::string name() {
//		return std::string("CWiseBinaryOp2<") + BinaryOp::name() + ", "
//				+ E1::name() + ", " + E2::name() + ", true>";
//	}
//#endif
//
//	/* ====================  MUTATORS      ======================================= */
//
//	/* ====================  OPERATORS     ======================================= */
//
//	STRONG_INLINE T operator[](size_t i) const {
//		return _op(_u[i], _v[i]);
//	}
//
//	STRONG_INLINE T operator()(size_t i, size_t j) const {
//		return _op(_u(i, j), _v(i, j));
//	}
//
//protected:
//	/* ====================  DATA MEMBERS  ======================================= */
//
//private:
//	/* ====================  DATA MEMBERS  ======================================= */
//	E1 const& _u;
//	Matrix<T2, major_order> _v;
//	Matrix<T, major_order> _tmpres;
//	bool _evaluted;
//
//	static const BinaryOp _op;
//
//};
/* ----------  end of template class CWiseBinaryOp2  ---------- */

//template<class BinaryOp, class E1, class E2, bool NeedTranspose =
//		traits<E1>::major_order != traits<E2>::major_order>
//struct BinaryOpReturnSelector {
//	typedef Matrix<typename traits<E2>::T, traits<E1>::major_order> return_type;
//	static CWiseBinaryOp2<BinaryOp, E1, return_type> getReturn(
//			MatExpression<E1> const& u, MatExpression<E2> const& v) {
//		return CWiseBinaryOp2<BinaryOp, E1, return_type>(u, v, true);
//	}
//};
//template<class BinaryOp, class E1, class E2>
//struct BinaryOpReturnSelector<BinaryOp, E1, E2, false> {
//	typedef E2 return_type;
//	static CWiseBinaryOp2<BinaryOp, E1, return_type> getReturn(
//			MatExpression<E1> const& u, MatExpression<E2> const& v) {
//		return CWiseBinaryOp2<BinaryOp, E1, return_type>(u, v);
//	}
//};
/* #define MAKE_CWISE_BINARY_OPERATOR(OPERATOR, CWISE_OP) template <class E1, class E2> \
 *             STRONG_INLINE CWiseBinaryOp2<CWISE_OP<typename traits<E1>::T,typename traits<E2>::T>, \
 *                    E1, typename BinaryOpReturnSelector<E1,E2>::return_type> \
 *             OPERATOR (MatExpression<E1> const& u, MatExpression<E2> const& v) { \
 *                 if(u.row() != v.row() || u.col() != v.col()) { \
 *                 }\
 *                 typedef typename BinaryOpReturnSelector<E1,E2>::return_type E3;\
 *                 typedef CWISE_OP<typename traits<E1>::T, typename traits<E2>::T> BinaryOp;\
 *                 if(traits<E1>::major_order == traits<E2>::major_order)\
 *                     return CWiseBinaryOp2<BinaryOp,E1,E3>(u,v);\
 *                 return CWiseBinaryOp2<BinaryOp,E1,E2>(u,E3(v));\
 *             }\
 */

template<class E1, class E2>
        STRONG_INLINE CWiseBinaryOp2<E1, E2> 

    operator+(MatExpression<E1> const& u, MatExpression<E2> const& v) {
        if (u.row() != v.row() || u.col() != v.col()) {
        }
//        typedef scalar_sum_op<typename traits<E1>::T, typename traits<E2>::T> BinaryOp;
        return CWiseBinaryOp2<E1, E2>(u, v);
    }

//        MAKE_CWISE_BINARY_OPERATOR(operator+,scalar_sum_op)
//        MAKE_CWISE_BINARY_OPERATOR(operator-,scalar_diff_op)
//        MAKE_CWISE_BINARY_OPERATOR(operator&,scalar_and_op)
//        MAKE_CWISE_BINARY_OPERATOR(operator|,scalar_or_op)
//        MAKE_CWISE_BINARY_OPERATOR(operator^,scalar_xor_op)

}
}

#endif
