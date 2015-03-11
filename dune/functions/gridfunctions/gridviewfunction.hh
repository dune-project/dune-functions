// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH

#include <memory>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>


namespace Dune {
namespace Functions {



template<class Signature, class GridView, template<class> class DerivativeTraits=DefaultDerivativeTraits, size_t bufferSize=64>
class GridViewFunction
{};



/**
 * \brief Wrapper class for functions defined on a GridView
 *
 * Being defined on a grid view means in particular that you can evaluate the function
 * in local coordinates of a given element of the grid view.
 *
 * \tparam GV The GridView that the function is defined on
 * \tparam Domain The domain type used for function arguments
 * \tparam Range The range type used for function values
 */
template<class Range, class Domain, class GV, template<class> class DerivativeTraits, size_t bufferSize>
class GridViewFunction<Range(Domain), GV, DerivativeTraits, bufferSize> :
  public GridFunction<Range(Domain), GridViewEntitySet<GV, 0>, DerivativeTraits, bufferSize>
{
  using Base = GridFunction<Range(Domain), GridViewEntitySet<GV, 0>, DerivativeTraits, bufferSize>;
public:
  using GridView = GV;

  using Base::Base;
};



} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH
