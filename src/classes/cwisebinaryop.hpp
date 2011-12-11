/*
 * =====================================================================================
 *
 *       Filename:  CWiseBinaryOp.hpp
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

#ifndef _CWISEBINARYOP_HPP_
#define _CWISEBINARYOP_HPP_

#ifdef __INTEL_COMPILER
    #include <boost/shared_ptr.hpp>
    using boost::shared_ptr;
#else
    #include <memory>
    using std::shared_ptr;
#endif

#include "matrix.hpp"
#include <utility>

//#define _op(A,B) (A)+(B)

namespace boss14420 {
namespace matrix {

template<class _BinaryOp, class _E1, class _E2>
struct traits<CWiseBinaryOp<_BinaryOp, _E1, _E2>> {
	typedef _BinaryOp BinaryOp;
	typedef _E1 E1;
	typedef _E2 E2;

	static const bool DifferenceMajorOrder = traits<E1>::major_order != traits<E2>::major_order;

        typedef const typename nested<E1,0>::type NestedE1;
        typedef const typename nested<E2,DifferenceMajorOrder ? 
            (NeedTranposeMask | NeedEvalute) : 0>::type NestedE2;

//        typedef const NestedE1 ConstNestedE1;
//        typedef const NestedE2 ConstNestedE2;

//        typedef const E1& NestedE1;
//        typedef const E2& NestedE2;

	typedef typename traits<E1>::T T1;
	typedef typename traits<E2>::T T2;

        typedef typename BinaryOp::T T;
	static const MajorOrder major_order = traits<E1>::major_order;

        typedef Matrix<T,major_order> ReturnType;
};

template<class BinaryOp, class E1, class E2, int Option>
    struct nested<CWiseBinaryOp<BinaryOp,E1,E2>, Option> {
        typedef CWiseBinaryOp<BinaryOp,E1,E2> Derived;
        typedef typename traits<Derived>::T T;
        static const MajorOrder major_order = traits<Derived>::major_order;

        typedef Matrix<typename traits<Derived>::T,
                        (Option & NeedTranposeMask) ? !major_order : major_order> type;
//        typedef CWiseBinaryOp<BinaryOp,E1,E2>& type;
    };

template<class BinaryOp, class E1, class E2>
    struct nested<CWiseBinaryOp<BinaryOp,E1,E2>, 0> {
        typedef CWiseBinaryOp<BinaryOp,E1,E2> type;
};

/*
 * =====================================================================================
 *        Class:  CWiseBinaryOp
 *  Description:  Matrix Add Expression class
 * =====================================================================================
 */
//template<class BinaryOp, class E1, class E2>
//class CWiseBinaryOp: public MatExpression<CWiseBinaryOp<BinaryOp, E1, E2>> {
//public:
//    	static const bool DifferenceMajorOrder = traits<E1>::major_order != traits<E2>::major_order;
//
//        typedef typename nested<E1,0>::type NestedE1;
//        typedef typename nested<E2,DifferenceMajorOrder ? 
//            NeedTranposeMask | NeedEvalute : 0>::type NestedE2;
//
////        typedef E1& NestedE1;
////        typedef E2& NestedE2;
//
//	typedef typename traits<E1>::T T1;
//	typedef typename traits<E2>::T T2;
//
//
//        typedef typename BinaryOp::T T;
//
//	static const MajorOrder major_order = traits<E1>::major_order;
//
//        typedef Matrix<T,major_order> ReturnType;
//
//       	/* ====================  LIFECYCLE     ======================================= */
//	CWiseBinaryOp(MatExpression<E1> const& u, MatExpression<E2> const& v) :
//			_u((NestedE1)u), _v((NestedE2)v)/*, _tmpres(), _evaluted(false) , _op(BinaryOp())*/ { /* constructor */
//#ifndef NDEBUG
//		std::cout << "Constructor CWiseBinaryOp ( " << MatExpression<E1>::name()
//				<< ", " << MatExpression<E2>::name() << " ) at " << this
//				<< " \n";
//		std::cout << "&u = " << &u << ", &v = " << &v << "\n";
//		std::cout << std::endl;
////                          std::cout << "u = " << u << "\nv = " << v << "\n\n";
//#endif
//
//	}
//
//	CWiseBinaryOp(const CWiseBinaryOp<BinaryOp, E1, E2>& mp) :
//			_u(mp._u), _v(mp._v)/* , _tmpres(mp._tmpres), 
//                        _evaluted(mp._evaluted) , _op( mp._op) */{
//#ifndef NDEBUG
//		std::cout << "Copy Constructor CWiseBinaryOp ( " << name() << " )\n\n";
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
//	ReturnType eval() const {
////		if (!_evaluted) {
//            ReturnType _tmpres(row(),col());
////			_tmpres.resize(row(), col());
//                        size_t i;
//#ifdef _OPENMP
//                        if(size() < MatExpression<CWiseBinaryOp<BinaryOp,E1,E2>>::ParallelThreshold)
//                            for(i = 0; i < size(); ++i) {
//				_tmpres[i] = _op(_u[i], _v[i]);
//                            }
//                        else
//#pragma omp parallel for private(i)
//#endif
//			for (i = 0; i < size(); ++i)
//				_tmpres[i] = _op(_u[i], _v[i]);
//
////		}
//#ifndef NDEBUG
//		std::cout << name() << " at " << this << std::endl;
//		std::cout << "_tmpres = " << _tmpres << std::endl << std::endl;
//#endif
//		return _tmpres;
//	}
//
//#ifndef NDEBUG
//	static std::string name() {
//		return std::string("CWiseBinaryOp<") + BinaryOp::name() + ", "
//				+ E1::name() + ", " + E2::name() + ">";
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
//        STRONG_INLINE operator ReturnType() const {
//            return eval();
//        }
//
//protected:
//	/* ====================  DATA MEMBERS  ======================================= */
//
//private:
//	/* ====================  DATA MEMBERS  ======================================= */
//	NestedE1 _u;
//	NestedE2 _v;
////	Matrix<T, major_order> _tmpres;
////	bool _evaluted;
//
//	static const BinaryOp _op;
//
//};
/* ----------  end of template class CWiseBinaryOp  ---------- */

/*
 * =====================================================================================
 *        Class:  CWiseBinaryOp
 *  Description:  Matrix Add Expression class
 * =====================================================================================
 */
template<class BinaryOp, class E1, class E2>
class CWiseBinaryOp : public MatExpression<CWiseBinaryOp<BinaryOp, E1, E2>> {
public:
	typedef typename traits<E1>::T T1;
	typedef typename traits<E2>::T T2;

//        typedef const E1& NestedE1;
//        typedef const E2& NestedE2;

	static const bool DifferenceMajorOrder = traits<E1>::major_order != traits<E2>::major_order;

        typedef typename nested<E1,0>::type NestedE1;
        typedef typename nested<E2,DifferenceMajorOrder ? 
            (NeedTranposeMask | NeedEvalute) : 0>::type NestedE2;

//        typedef NestedE1 ConstNestedE1;
//        typedef NestedE2 ConstNestedE2;

	static const MajorOrder major_order = traits<E1>::major_order;
	typedef decltype(BinaryOp()(T1(),T2())) T;

        typedef Matrix<T,major_order> ReturnType;

	/* ====================  LIFECYCLE     ======================================= */

	CWiseBinaryOp(MatExpression<E1> const& u, MatExpression<E2> const& v)
//            : _u(u),
//              _v(v)/*,
            : _u(pass_to_expression<NestedE1>()(u)),
              _v(pass_to_expression<NestedE2>()(v))
               /*,_tmpres(), _evaluted(false) */{
	    		
#ifndef NDEBUG
		std::cout << "Constructor CWiseBinaryOp ( " << MatExpression<E1>::name()
				<< ", " << MatExpression<E2>::name() << ") at " << this
				<< " \n";
		std::cout << "&u = " << &u << ", &v = " << &v << "\n";
		std::cout << "&_u = " << &_u << ", &_v = " << &_v << "\n";
		std::cout << std::endl;
		//                          std::cout << "u = " << u << "\nv = " << v << "\n\n";
#endif
	}

	CWiseBinaryOp(const CWiseBinaryOp<BinaryOp, E1, E2>& mp) :
			_u(mp._u), _v(mp._v)/* , _tmpres(mp._tmpres), 
                        _evaluted(mp._evaluted)*//*, _op( mp._op) */{
#ifndef NDEBUG
		std::cout << "Copy Constructor CWiseBinaryOp ( " << name() << " )\n\n";
#endif
	}

        CWiseBinaryOp(CWiseBinaryOp<BinaryOp, E1, E2>&& mp) :
                _u(pass_to_expression<NestedE1>()(mp._u)), 
                _v(pass_to_expression<NestedE2>()(mp._v))/*, 
                _tmpres(std::move(mp._tmpres)), _evaluted(mp._evaluted)*/{
#ifndef NDEBUG
                    std::cout << "Move Constructor CWiseBinaryOp ( " << name() << " )\n\n";
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

//	const Matrix<T, major_order>& eval() {
//        ReturnType eval() {
//		if (!_evaluted) {
////                    _ptrres.reset(new ReturnType(row(),col()));
//			_tmpres.resize(row(), col());
////                    ReturnType& _tmpres = *_ptrres;
//                        size_t i;
//#ifdef _OPENMP
//                        if(size() < MatExpression<CWiseBinaryOp<BinaryOp,E1,E2>>::ParallelThreshold)
//                            for(i = 0; i < size(); ++i) {
//				_tmpres[i] = _op(_u[i], _v[i]);
//                            }
//                        else
//#pragma omp parallel for private(i)
//#endif
//			for (i = 0; i < size(); ++i)
//				_tmpres[i] = _op(_u[i], _v[i]);
//                    _evaluted = true;
//		}
//#ifndef NDEBUG
//		std::cout << name() << " at " << this << std::endl;
//		std::cout << "_tmpres = " << _tmpres << std::endl << std::endl;
//#endif
//		return _tmpres;
//	}

#ifndef NDEBUG
	static std::string name() {
		return std::string("CWiseBinaryOp<") + BinaryOp::name() + ", "
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

//        operator ReturnType () {
//            return eval();
//        }

protected:
	/* ====================  DATA MEMBERS  ======================================= */

private:
	/* ====================  DATA MEMBERS  ======================================= */
	NestedE1 _u;
//	Matrix<T2, major_order> _v;
        NestedE2 _v;
//	shared_ptr<ReturnType> _ptrres;
//        ReturnType _tmpres;
//	bool _evaluted;

	static const BinaryOp _op;

};
/* ----------  end of template class CWiseBinaryOp  ---------- */

//template<class BinaryOp, class E1, class E2, bool NeedTranspose =
//		traits<E1>::major_order != traits<E2>::major_order>
//struct BinaryOpReturnSelector {
//	typedef Matrix<typename traits<E2>::T, traits<E1>::major_order> return_type;
//	static CWiseBinaryOp<BinaryOp, E1, return_type> getReturn(
//			MatExpression<E1> const& u, MatExpression<E2> const& v) {
//		return CWiseBinaryOp<BinaryOp, E1, return_type>(u, v, true);
//	}
//};
//template<class BinaryOp, class E1, class E2>
//struct BinaryOpReturnSelector<BinaryOp, E1, E2, false> {
//	typedef E2 return_type;
//	static CWiseBinaryOp<BinaryOp, E1, return_type> getReturn(
//			MatExpression<E1> const& u, MatExpression<E2> const& v) {
//		return CWiseBinaryOp<BinaryOp, E1, return_type>(u, v);
//	}
//};
/* #define MAKE_CWISE_BINARY_OPERATOR(OPERATOR, CWISE_OP) template <class E1, class E2> \
 *             STRONG_INLINE CWiseBinaryOp<CWISE_OP<typename traits<E1>::T,typename traits<E2>::T>, \
 *                    E1, typename BinaryOpReturnSelector<E1,E2>::return_type> \
 *             OPERATOR (MatExpression<E1> const& u, MatExpression<E2> const& v) { \
 *                 if(u.row() != v.row() || u.col() != v.col()) { \
 *                 }\
 *                 typedef typename BinaryOpReturnSelector<E1,E2>::return_type E3;\
 *                 typedef CWISE_OP<typename traits<E1>::T, typename traits<E2>::T> BinaryOp;\
 *                 if(traits<E1>::major_order == traits<E2>::major_order)\
 *                     return CWiseBinaryOp<BinaryOp,E1,E3>(u,v);\
 *                 return CWiseBinaryOp<BinaryOp,E1,E2>(u,E3(v));\
 *             }\
 */

template<class E1, class E2>
        STRONG_INLINE CWiseBinaryOp<scalar_sum_op<typename traits<E1>::T, 
            typename traits<E2>::T>, E1, E2> 

    operator+(MatExpression<E1> const& u, MatExpression<E2> const& v) {
        if (u.row() != v.row() || u.col() != v.col()) {
        }
        typedef scalar_sum_op<typename traits<E1>::T, typename traits<E2>::T> BinaryOp;
        return CWiseBinaryOp<BinaryOp, E1, E2>(u, v);
    }

//        MAKE_CWISE_BINARY_OPERATOR(operator+,scalar_sum_op)
//        MAKE_CWISE_BINARY_OPERATOR(operator-,scalar_diff_op)
//        MAKE_CWISE_BINARY_OPERATOR(operator&,scalar_and_op)
//        MAKE_CWISE_BINARY_OPERATOR(operator|,scalar_or_op)
//        MAKE_CWISE_BINARY_OPERATOR(operator^,scalar_xor_op)

}
}

#endif
