// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_LOCALDERIVATIVE_TRAITS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_LOCALDERIVATIVE_TRAITS_HH


#include <dune/functions/common/defaultderivativetraits.hh>


namespace Dune {
namespace Functions {


/**
 * \brief Derivative traits for local functions
 *
 * \ingroup FunctionUtility
 *
 * This provides derivative traits to be used by local functions with
 * given EntitySet and global DerivativeTraits.
 *
 * The reason why this is needed is that local functions return
 * derivatives wrt to global coordinates although the arguments
 * for these derivatives are local coordinates.
 */
template<class EntitySet, template<class> class DerivativeTraits=DefaultDerivativeTraits>
struct LocalDerivativeTraits
{
  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Domain = typename EntitySet::GlobalCoordinate;

  template<class Signature>
  struct Traits
  {
    typedef InvalidRange Range;
  };

  template<class R>
  struct Traits<R(LocalDomain)>
  {
    using Range = typename DerivativeTraits<R(Domain)>::Range;
  };
};


}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_LOCALDERIVATIVE_TRAITS_HH
