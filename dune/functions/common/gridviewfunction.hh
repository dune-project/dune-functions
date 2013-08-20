// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH

#include <memory>
#include <dune/functions/common/differentiablefunction.hh>

namespace Dune {

namespace Functions {



template<typename GV, typename RT>
class GridViewFunction
  : public DifferentiableFunction<typename GV::template Codim<0>::Geometry::GlobalCoordinate,RT>
{

  typedef DifferentiableFunction<
    typename GV::template Codim<0>::Geometry::GlobalCoordinate,
    RT
    > Base;

  typedef GV GridViewType;
  typedef typename Base::DomainType DomainType;
  typedef typename Base::RangeType RangeType;
  typedef typename Base::DerivativeRangeType DerivativeRangeType;

  typedef GridViewFunction<GV,DerivativeRangeType> DerivativeType;


  typedef typename GV::template Codim<0>::Geometry::LocalCoordinate LocalDomainType;
  typedef typename GV::template Codim<0>::Entity ElementType;

  virtual void evaluate(const ElementType& e, const LocalDomainType& coord, RangeType& r) const = 0;

  virtual DerivativeType* derivative() const = 0;

};



} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH
