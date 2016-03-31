// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_FUNCTION_FROM_CALLABLE_HH
#define DUNE_FUNCTIONS_COMMON_FUNCTION_FROM_CALLABLE_HH

#include <dune/common/function.hh>

#include <dune/functions/common/signature.hh>


namespace Dune {
namespace Functions {



template<class Signature, class F,
  class FunctionInterface = typename Dune::VirtualFunction<
    typename SignatureTraits<Signature>::RawDomain,
    typename SignatureTraits<Signature>::RawRange> >
class FunctionFromCallable;

/**
 * \brief Wrap a callable object as Dune::Function or Dune::VirtualFunction
 *
 * \ingroup FunctionUtility
 *
 * You can use this to implement a DifferentiableFunction including
 * a variable number of derivatives using callable objects. All
 * types that can be assigned to std::function are supported,
 * i.e. functions, functors, lambdas, ...
 *
 * \tparam Signature Function signature
 * \tparam F Type of wrapped function
 * \tparam FunctionInterface Interface to implement, this can be either Dune::Function or Dune::VirtualFunction
 */
template<class Range, class Domain, class F, class FunctionInterface>
class FunctionFromCallable<Range(Domain), F, FunctionInterface> :
  public FunctionInterface
{
public:

  /**
   * \brief Create VirtualFunction from callable object
   *
   * This will store a std::function<Domain,Range> wrapping the callable
   * and pass evaluate() to this object.
   *
   * \tparam F Anything that can be assigned to std::function<Domain, Range> (funcions, functors, lambdas,...)
   * \param f Callable object to use for evaluate()
   */
  FunctionFromCallable(F&& f) :
    f_(f)
  {}

  FunctionFromCallable(const F& f) :
    f_(f)
  {}

  /**
   * \brief Evaluate function
   *
   * This call is passed to the function
   */
  void evaluate(const Domain& x, Range&y) const
  {
    y = f_(x);
  }

private:
  F f_;
};



} // namespace Functions
} // namespace Dune

#endif //DUNE_FUNCTIONS_COMMON_FUNCTION_FROM_CALLABLE_HH
