// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_FUNCTION_FROM_CALLABLE_HH
#define DUNE_FUNCTIONS_COMMON_FUNCTION_FROM_CALLABLE_HH


#include <memory>
#include <functional>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
#include <dune/common/std/final.hh>
#else
 #ifndef DUNE_FINAL
  #define DUNE_FINAL
 #endif
#endif
#include <dune/common/exceptions.hh>
#include <dune/functions/common/differentiablefunction.hh>


namespace Dune {
namespace Functions {


template<class D, class R>
class FunctionFromCallable;

// Helper class that exports the type for FunctionFromCallable<D,R>::Derivative
// and constructs the corresponding object if possible.
// Notice that the second template argument is the DerivativeRange and NOT the Range
template<class D, class DR>
struct FunctionFromCallableDerivativeTraits
{
  typedef FunctionFromCallable<D, DR> Derivative;

  template<class... DF>
  static std::shared_ptr<Derivative> makeIfImplemented(const DF&... df)
  {
    return std::make_shared<Derivative>(df...);
  }
};

// Specialization for DerivativeRange==InvalidRange. In this case
// FunctionFromCallable<D,R>::Derivative is the interface class.
template<class D>
struct FunctionFromCallableDerivativeTraits<D, InvalidRange>
{
  typedef DifferentiableFunction<D, InvalidRange> Derivative;

  template<class... DF>
  static std::shared_ptr<Derivative> makeIfImplemented(const DF&... df)
  {
    DUNE_THROW(Dune::NotImplemented, "DerivativeRange is not implemented for this type.");
  }
};



/**
 * \brief Wrap a callable object as Dune::DifferentiableFunction
 *
 * You can use this to implement a DifferentiableFunction including
 * a variable number of derivatives using callable objects. All
 * types that can be assigned to std::function are supported,
 * i.e. functions, functors, lambdas, ...
 *
 * \tparam D The domain type
 * \tparam R The range type
 *
 */
template<class D, class R>
class FunctionFromCallable :
  public DifferentiableFunction<D, R>
{
    typedef DifferentiableFunction<D, R> Base;

  public:

    typedef typename Base::Range Range;
    typedef typename Base::Domain Domain;
    typedef typename Base::DerivativeRange DerivativeRange;

    typedef typename FunctionFromCallableDerivativeTraits<Domain, DerivativeRange>::Derivative Derivative;

    /**
     * \brief Create function from callable object
     *
     * This will store a std::function<Domain,Range> wrapping the callable
     * and pass evaluate() to this object.
     *
     * \tparam F Anything that can be assigned to std::function<Domain, Range> (funcions, functors, lambdas,...)
     * \param f Callable object to use for evaluate()
     */
    template<typename F>
    FunctionFromCallable(const F& f)
    {
      f_ = f;
    }

    /**
     * \brief Create function and derivatives from callable objects
     *
     * This will store a std::function<Domain,Range> wrapping the callable
     * and pass evaluate() to this object.
     *
     * All further arguments are used to implement derivatives.
     * The n-th argument should be assignable to std::function<Domain, NDR>
     * with NDR being the Range of the n-th derivative.
     * By passing n arguments the created function will support
     * the evaluation of n derivatives.
     *
     * \tparam F Anything that can be assigned to std::function<Domain, Range> (funcions, functors, lambdas,...)
     * \tparam DF Variable length arg list for derivatives
     *
     * \param f Callable object to use for evaluate()
     * \param df Variable length arg list of Callable objects to use for evaluate() of derivatives
     */
    template<typename F, typename... DF>
    FunctionFromCallable(const F& f, const DF&... df)
    {
      f_ = f;
      derivative_ = FunctionFromCallableDerivativeTraits<Domain, DerivativeRange>::makeIfImplemented(df...);
    }

    /**
     * \brief Evaluate function
     *
     * This call is passed to the stored std::function
     */
    virtual void evaluate(const Domain& x, Range&y) const DUNE_FINAL
    {
      y = f_(x);
    }

    /**
     * \brief Return derivative
     *
     * If a proper number of derivatives where provided
     * this will return the derivative as FunctionFromCallable<Domain, DerivativeRange>
     */
    virtual std::shared_ptr<typename Base::Derivative> derivative() const DUNE_FINAL
    {
      if (not derivative_)
        DUNE_THROW(Dune::NotImplemented, "This derivative was not provided.");
      return derivative_;
    }

  private:
    std::function<R(D)> f_;
    std::shared_ptr<typename Base::Derivative> derivative_;
};




} // namespace Functions
} // namespace Dune

#endif //DUNE_FUNCTIONS_COMMON_FUNCTION_FROM_CALLABLE_HH
