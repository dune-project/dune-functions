// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
#include <memory>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/function.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/final.hh>

namespace Dune {

namespace Functions {



template<class DT, class RT> class DifferentiableFunction;

template<class FImp> class FunctionHandle;



class InvalidRange
{};


template<class DT, class RT>
struct DerivativeTraits
{
  typedef InvalidRange DerivativeRange;
};

template<>
struct DerivativeTraits<double, double>
{
  typedef double DerivativeRange;
};

template<typename F, int n>
struct DerivativeTraits<FieldVector<F,n>,F>
{
  typedef FieldVector<F,n> DerivativeRange;
};

template<typename F, int n, int m>
struct DerivativeTraits<FieldVector<F,n>,FieldVector<F,m> >
{
  typedef FieldMatrix<F,n,m> DerivativeRange;
};

template<typename F, int n, int m>
struct DerivativeTraits<FieldVector<F,n>,FieldMatrix<F,m,1> >
{
  typedef FieldMatrix<F,n,m> DerivativeRange;
};



/** \brief Abstract base class for functions that allow to compute a derivative
 *
 * \tparam DT Type used for the independent variable (may be a vector type for functions of several variables)
 * \tparam RT Type used for function values
 *
 * We view differentiation as a map between functions.  In other words, the derivative of a function
 * is again a function.  The derivative function implements the same DifferentiableFunction base class,
 * only the type for function values will have changed.
 */
template<class DT, class RT>
class DifferentiableFunction :
    public Dune::VirtualFunction<DT, RT>
{
  public:
    /** \brief Type used for the independent variable (may be a vector type for functions of several variables) */
    typedef DT Domain;

    /** \brief Type used for function values */
    typedef RT Range;

    /** \brief Type of the values of the derivative function */
    typedef typename DerivativeTraits<DT, RT>::DerivativeRange DerivativeRange;

    /** \brief Type of the derivative function */
    typedef DifferentiableFunction<Domain, DerivativeRange> Derivative;

    /** \brief Get the derivative function
     *
     * ### Life time of the return value
     *
     * This method returns a std::shared_ptr<Derivative>. Since covariant return
     * values do not work with shared_ptrs the type 'Derivative' will also be the
     * interface class for actual implementations. In order to avoid virtual method
     * calls if the implementation type is known one should always use the free
     * method derivative(f) instead of f.derivative(). This will cast the shared_ptr
     * to the appropriate implementation class.
     */
    virtual std::shared_ptr<Derivative> derivative() const = 0;
};


/** \brief Get the derivative of a DifferentiableFunction
 *
 * \tparam FImp A DifferentiableFunction, i.e., a function that has a derivative() method
 *
 * Using this methods avoids ever having to handle the C-style pointer that is returned
 * by DifferentiableFunction::derivative.
 */
template<class FImp>
std::shared_ptr<typename FImp::Derivative> derivative(const FImp& f)
{
  return std::static_pointer_cast<typename FImp::Derivative>(f.derivative());
}

template<class FImp>
std::shared_ptr<typename FImp::Derivative> derivative(const std::shared_ptr<FImp>& f)
{
  return std::static_pointer_cast<typename FImp::Derivative>(f->derivative());
}

} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
