// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH

#include <type_traits>

#include <dune/functions/common/signature.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>


namespace Dune {
namespace Functions {

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
  using Geometry = typename Element::Geometry;

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
//    geometry_ = element_.geometry();
  }

  void unbind()
  {}

  Range operator()(const LocalDomain& x) const
  {
//    return f_(geometry_.global(x));
    return f_(element_.geometry().global(x));
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

//  Geometry geometry_;
  Element element_;
  F f_;
};



template<class Signature, class GV, class F, template<class> class DerivativeTraits=DefaultDerivativeTraits>
class AnalyticGridViewFunction;


/**
 * \brief Class wrapping any differentiable function as grid function
 *
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
  using LocalFunction = LocalAnalyticGridViewFunction<Range(LocalDomain), GridView, F, DerivativeTraits>;

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



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH
