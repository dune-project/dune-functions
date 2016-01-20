// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH

#include <type_traits>

#include <dune/common/typeutilities.hh>

#include <dune/functions/common/signature.hh>
#include <dune/functions/common/optional.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/localderivativetraits.hh>


namespace Dune {
namespace Functions {

namespace Imp {

template<class Signature, class GV, class FLocal, template<class> class DerivativeTraits=DefaultDerivativeTraits>
class LocalAnalyticGridViewFunction;

template<class Range, class LocalDomain, class GV, class F, template<class> class DerivativeTraits>
class LocalAnalyticGridViewFunction<Range(LocalDomain), GV, F, DerivativeTraits>
{
public:
  using Signature = Range(LocalDomain);
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(LocalDomain);

  using GridView = GV;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;
//  using Geometry = typename Element::Geometry;
  using Geometry = typename std::decay<typename Element::Geometry>::type;

  // Use the inderiction via derivativeIfImplemented to also support
  // function types F that do not implement derivative. In this case
  // the interface type DifferentiableFunction is used a dummy for
  // the derivative type
  using DerivativeDummy = DifferentiableFunction<DerivativeSignature>;
  using GlobalRawDerivative = decltype(Imp::derivativeIfImplemented<DerivativeDummy, F>(std::declval<F>()));
  using LocalDerivative = LocalAnalyticGridViewFunction<DerivativeSignature, GridView, GlobalRawDerivative, DerivativeTraits>;

  template<class FT, disableCopyMove<LocalAnalyticGridViewFunction, FT> = 0>
  LocalAnalyticGridViewFunction(FT&& f) :
    f_(std::forward<FT>(f))
  {}


  void bind(const Element& element)
  {
    element_ = element;

    // We'd like to do
    //
    //   geometry_ = element_.geometry();
    //
    // But since Geometry is not assignable we
    // have to reconstruct it - argh
    geometry_.release();
    geometry_.emplace(element_.geometry());
  }

  void unbind()
  {}

  Range operator()(const LocalDomain& x) const
  {
    return f_(geometry_.value().global(x));
  }

  const Element& localContext() const
  {
    return element_;
  }

  friend LocalDerivative derivative(const LocalAnalyticGridViewFunction& t)
  {
    return LocalDerivative(Imp::derivativeIfImplemented<DerivativeDummy, F>(t.f_));
  }

private:

  // Hack around the fact that Geometry is not default constructible.
  Optional<Geometry> geometry_;

  Element element_;
  F f_;
};

} // end namespace Imp




template<class Signature, class GV, class F, template<class> class DerivativeTraits=DefaultDerivativeTraits>
class AnalyticGridViewFunction;


/**
 * \brief Class wrapping any differentiable function as grid function
 *
 * \ingroup FunctionImplementations
 */
template<class Range, class Domain, class GV, class F, template<class> class DerivativeTraits>
class AnalyticGridViewFunction<Range(Domain), GV, F, DerivativeTraits>
{
public:
  using Signature = Range(Domain);
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

  using GridView = GV;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;
  using Geometry = typename Element::Geometry;

  // Use the inderiction via derivativeIfImplemented to also support
  // function types F that do not implement derivative. In this case
  // the interface type DifferentiableFunction is used a dummy for
  // the derivative type
  using DerivativeDummy = DifferentiableFunction<DerivativeSignature>;
  using GlobalRawDerivative = decltype(Imp::derivativeIfImplemented<DerivativeDummy, F>(std::declval<F>()));
  using Derivative = AnalyticGridViewFunction<DerivativeSignature, GridView, GlobalRawDerivative, DerivativeTraits>;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using LocalFunction = typename Imp::LocalAnalyticGridViewFunction<Range(LocalDomain), GridView, F, LocalDerivativeTraits<EntitySet, DerivativeTraits>::template Traits>;

  template<class FT>
  AnalyticGridViewFunction(FT&& f, const GridView& gridView) :
    f_(std::forward<FT>(f)),
    entitySet_(gridView)
  {}

  Range operator()(const Domain& x) const
  {
    return f_(x);
  }

  friend Derivative derivative(const AnalyticGridViewFunction& t)
  {
    return Derivative(Imp::derivativeIfImplemented<DerivativeDummy, F>(t.f_), t.entitySet_.gridView());
  }

  friend LocalFunction localFunction(const AnalyticGridViewFunction& t)
  {
    return LocalFunction(t.f_);
  }

  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:
  F f_;
  EntitySet entitySet_;
};



/**
 * \brief Construct AnalyticGridViewFunction from function and grid view
 *
 * The returned function supports localFunction() and stores a copy of the
 * original function.
 *
 * \param f A function object supporting argument compatible with global coordinates
 * \param gridView The GridView the function should act on.
 *
 * \returns A function that models the GridFunction interface.
 */
template<class F, class GridView>
AnalyticGridViewFunction<
  typename std::result_of<F(typename GridView::template Codim<0>::Geometry::GlobalCoordinate)>::type  // Range
  (typename GridView::template Codim<0>::Geometry::GlobalCoordinate),                                 // Domain
  GridView,
  typename std::decay<F>::type >                                                                      // Raw type of F (without & or &&)
  makeAnalyticGridViewFunction(F&& f, const GridView& gridView)
{
  using Domain = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
  using Range = typename std::result_of<F(Domain)>::type;
  using FRaw = typename std::decay<F>::type;

  return AnalyticGridViewFunction<Range(Domain), GridView, FRaw>(std::forward<F>(f), gridView);
}



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH
