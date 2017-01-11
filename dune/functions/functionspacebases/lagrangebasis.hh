// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH

#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>


namespace Dune {
namespace Functions {



namespace BasisBuilder {

template<std::size_t k>
Imp::PQkNodeFactoryBuilder<k> lagrange()
{
  return{};
}

} // end namespace BasisBuilder



/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of PQkNodeFactory.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using LagrangeBasis = DefaultGlobalBasis<PQkNodeFactory<GV, k, FlatMultiIndex<std::size_t>> >;





} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
