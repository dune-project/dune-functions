// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_COMPOSEDGRIDFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_COMPOSEDGRIDFUNCTION_HH

#include <type_traits>
#include <tuple>

#include <dune/common/typeutilities.hh>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/referencehelper.hh>


namespace Dune {
namespace Functions {



/**
 * \brief Composition of grid functions with another function
 *
 * \ingroup FunctionImplementations
 *
 * For given inner grid functions g0, ..., gn and an
 * outer function f this creates a grid function
 * representing f(g0(x), ..., gn(x)). The only assumption
 * made, is that the range types of the inner functions
 * can be passed to the outer ones, and that all grid
 * functions are defined on the same EntitySet.
 *
 * Notice that all functions are captured by value.
 * To store references you can pass use std::ref().
 *
 * \tparam OF Type of outer function. std::reference_wrapper is supported.
 * \tparam IF Types of inner outer functions. std::reference_wrapper is supported.
 */
template<class OF, class... IF>
class ComposedGridFunction
{
  using InnerFunctions = std::tuple<IF...>;
  using InnerLocalFunctions = std::tuple<decltype(localFunction(resolveRef(std::declval<const IF&>())))...>;

  template<std::size_t i>
  using InnerFunction = std::decay_t<ResolveRef_t<std::tuple_element_t<i, InnerFunctions>>>;

  using OuterFunction = OF;

public:

  using EntitySet = typename InnerFunction<0>::EntitySet;
  using Element = typename EntitySet::Element;

  using Domain = typename EntitySet::GlobalCoordinate;
  using LocalDomain = typename EntitySet::LocalCoordinate;

  using Range = decltype(std::declval<OF>()(std::declval<IF>()(std::declval<Domain>())...));

private:

  using Traits = Imp::GridFunctionTraits<Range(Domain), EntitySet, DefaultDerivativeTraits, 16>;

  class LocalFunction
  {
  public:

    LocalFunction(const ComposedGridFunction& globalFunction) :
      globalFunction_(globalFunction),
      innerLocalFunctions_(globalFunction.innerLocalFunctions())
    {}


    void bind(const Element& element)
    {
      std::apply([&](auto&... innerFunction) {
          (innerFunction.bind(element),...);
      }, innerLocalFunctions_);
    }

    void unbind()
    {
      std::apply([&](auto&... innerFunction) {
          (innerFunction.unbind(),...);
      }, innerLocalFunctions_);
    }

    Range operator()(const LocalDomain& x) const
    {
      return std::apply([&](const auto&... innerFunction) {
          return globalFunction_.outerFunction_(innerFunction(x)...);
      }, innerLocalFunctions_);
    }

    const Element& localContext() const
    {
      return std::get<0>(innerLocalFunctions_).localContext();
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:
    const ComposedGridFunction& globalFunction_;
    InnerLocalFunctions innerLocalFunctions_;
  };

public:

  /**
   * \brief Create ComposedGridFunction
   *
   *
   * Outer and inner functions will be captured by value.
   * To store references you can pass use std::ref().
   *
   * \param outerFunction The outer function to be composed with the grid functions.
   * \param innerFunctions The inner grid functions
   */
  template<class OFT, class... IFT,
    disableCopyMove<ComposedGridFunction, OFT> = 0>
  ComposedGridFunction(OFT&& outerFunction, IFT&&... innerFunctions) :
    outerFunction_(std::forward<OFT>(outerFunction)),
    innerFunctions_(std::forward<IFT>(innerFunctions)...)
  {}

  Range operator()(const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend typename Traits::DerivativeInterface derivative(const ComposedGridFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend LocalFunction localFunction(const ComposedGridFunction& cgf)
  {
    return LocalFunction(cgf);
  }

  const EntitySet& entitySet() const
  {
    return resolveRef(std::get<0>(innerFunctions_)).entitySet();
  }

protected:

  InnerLocalFunctions innerLocalFunctions() const
  {
    return std::apply([&](const auto&... innerFunction) {
        return std::make_tuple(localFunction(resolveRef(innerFunction))...);
    }, innerFunctions_);
  }

  OuterFunction outerFunction_;
  InnerFunctions innerFunctions_;
};



/**
 * \brief Compose grid functions with another function
 *
 * \ingroup FunctionImplementations
 *
 * For given inner grid functions g0, ..., gn and an
 * outer function f this creates a grid function
 * representing f(g0(x), ..., gn(x)). The only assumption
 * made, is that the range types of the inner functions
 * can be passed to the outer ones, and that all grid
 * functions are defined on the same EntitySet.
 *
 * Notice that all functions are captured by value.
 * To store references you can pass use std::ref().
 *
 * \param outerFunction The outer function to be composed with the grid functions.
 * \param innerFunctions The inner grid functions
 *
 * \returns A grid function defined on the same EntitySet as the input functions.
 */
template<class OF, class... IF>
auto makeComposedGridFunction(OF&& outerFunction, IF&&... innerFunction)
{
  using ComposedGridFunctionType = ComposedGridFunction<std::decay_t<OF>, std::decay_t<IF>...>;
  return ComposedGridFunctionType(std::forward<OF>(outerFunction), std::forward<IF>(innerFunction)...);
}



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_COMPOSEDGRIDFUNCTION_HH
