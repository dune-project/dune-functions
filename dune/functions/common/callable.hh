// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_CALLABLE_HH
#define DUNE_FUNCTIONS_COMMON_CALLABLE_HH

#include <memory>
#include <functional>

#include <dune/common/function.hh>
#include <dune/common/shared_ptr.hh>


namespace Dune {
namespace Functions {




/**
 * \brief Wrap a Dune::VirtualFunction into a callable object
 *
 * \ingroup FunctionImplementations
 *
 * \tparam F Some function deriving from Dune::VirtualFunction
 *
 * This class provides an operator() that forwards to the evaluate
 * method. In order to allow this F::RangeType must be default constructible.
 * This class is copyable and ownership policy of the wrapped function
 * (see constructors) will be the same for all copies.
 *
 * This models the \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> concept.
 *
 * Using this wrapper you can e.g. pass a Dune::VirtualFunction<D,R> to
 * the standard wrapper std::function<R(D)>.
 */
template<class F>
class CallableFunctionWrapper
{
    typedef typename F::Range Range;
    typedef typename F::Domain Domain;

  public:

    /**
     * \brief Instanciate from reference to f
     *
     * The CallableFunctionWrapper will not take ownership
     * of the provided function.
     */
    CallableFunctionWrapper(const F& f)
    {
      f_ = Dune::stackobject_to_shared_ptr(f);
    }

    /**
     * \brief Instanciate from shared_ptr to f
     *
     * The CallableFunctionWrapper will share ownership
     * with the provided function.
     */
    CallableFunctionWrapper(const std::shared_ptr<const F>& f) :
      f_(f)
    {}

    /**
     * \brief Forward operator() to F::evaluate()
     *
     * This uses the default constructor of F::RangeType
     */
    Range operator()(const Domain& x) const
    {
      Range y;
      f_->evaluate(x, y);
      return y;
    }

  private:
    std::shared_ptr<const F> f_;
};



/**
 * \brief Create a callable object from some Dune::VirtualFunction
 *
 * \ingroup FunctionImplementations
 *
 * \tparam F Function type derived from Dune::VirtualFunction
 * \param f The function to be wrapper
 *
 * The returned object will only be valid as long f is valid.
 * You can e.g. do the following:
 * \code
 *   // Create some F derived from VirtualFunction<D,R>
 *   F f;
 *
 *   // store callable directly
 *   auto fc = callable(f);
 *
 *   // store callable through default wrapper
 *   std::function<R(D)> stdF = callable(f)
 * \endcode
 */
template<class F>
CallableFunctionWrapper<F> callable(const F& f)
{
  return CallableFunctionWrapper<F>(f);
}


/**
 * \brief Create a callable object from shared_ptr<F>
 *
 * \ingroup FunctionImplementations
 *
 * \tparam F Function type derived from Dune::VirtualFunction
 * \param fp shared_ptr<F> to the function to be wrapper
 *
 * The returned object will share ownership of fp
 * using a shared_ptr. You can e.g. do the following:
 * \code
 *   // Create some F derived from VirtualFunction<D,R>
 *   auto f = make_shared<F>();
 *
 *   // store callable directly
 *   auto f1 = callable(f);
 *
 *   // store callable through default wrapper
 *   std::function<R(D)> f2 = callable(f)
 *
 *   // Create some F derived from VirtualFunction<D,R> and only store
 *   // it in the callable wrapper
 *   auto f3 = callable(make_shared<F>());
 * \endcode
 */
template<class F>
CallableFunctionWrapper<F> callable(const std::shared_ptr<F>& fp)
{
  return CallableFunctionWrapper<F>(fp);
}



} // namespace Functions
} // namespace Dune

#endif //DUNE_FUNCTIONS_COMMON_CALLABLE_HH
