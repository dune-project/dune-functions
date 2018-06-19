// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>


namespace Dune {
namespace Functions {



/**
 * \brief A pre-basis for PQ-lagrange bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions
 * \tparam MI  Type to be used for multi-indices
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 */
template<typename GV, int k, class MI>
using PQkPreBasis
  DUNE_DEPRECATED_MSG("PQkPreBasis is deprecated. Use LagrangePreBasis instead.")
  = LagrangePreBasis<GV, k, MI>;

template<typename GV, int k, typename TP>
using PQkNode
  DUNE_DEPRECATED_MSG("PQkNode is deprecated. Use LagrangeNode instead.")
  = LagrangeNode<GV, k, TP>;



namespace BasisFactory {



/**
 * \brief Create a pre-basis factory that can create a PQ_k pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of ansatz functions
 */
template<std::size_t k>
auto DUNE_DEPRECATED_MSG("BasisFactory::pq is deprecated. Use BasisFactory::lagrange instead.")
  pq()
{
  return Imp::LagrangePreBasisFactory<k>();
}

} // end namespace BasisFactory


// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

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
 * of PQkPreBasis.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using PQkNodalBasis
  DUNE_DEPRECATED_MSG("PQkNodalBasis is deprecated. Use LagrangeBasis instead.")
  = DefaultGlobalBasis<LagrangePreBasis<GV, k, FlatMultiIndex<std::size_t>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
