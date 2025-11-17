// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>




namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   LagrangeDGPreBasis
//   LagrangeDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename R=double>
using LagrangeDGNode = LagrangeNode<GV, k, R>;

/** \brief PreBasis implementation for a Lagrangean-DG finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The order of ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   Range type used for shape function values
 */
template<typename GV, int k, typename R=double>
class LagrangeDGPreBasis :
  public Dune::Functions::LeafPreBasisMapperMixin<GV>
{
  using Base = Dune::Functions::LeafPreBasisMapperMixin<GV>;

  static constexpr bool useDynamicOrder = (k<0);

  static MCMGLayout dofLayout(int order)
  {
    return [order](Dune::GeometryType type, size_t dimGrid) {
      if (type.dim() == dimGrid)
      {
        if (type.isLine())
          return order+1;
        if (type.isTriangle())
          return (order+1)*(order+2)/2;
        if (type.isQuadrilateral())
          return (order+1)*(order+1);
        if (type.isTetrahedron())
          return (order+1)*(order+2)*(order+3)/6;
        if (type.isPrism())
          return (order+1)*(order+1)*(order+2)/2;
        if (type.isHexahedron())
          return (order+1)*(order+1)*(order+1);
        if (type.isPyramid())
          return (order+1)*(order+2)*(2*order+3)/6;
        DUNE_THROW(Dune::NotImplemented, "Using LagrangeDGPreBasis with non-supported GeometryType");
      }
      return 0;
    };
  }

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = typename Base::size_type;

  // Precompute the number of dofs per entity type
  [[deprecated("This constant will be removed after Dune 2.11")]] const static int dofsPerEdge        = k+1;
  [[deprecated("This constant will be removed after Dune 2.11")]] const static int dofsPerTriangle    = (k+1)*(k+2)/2;
  [[deprecated("This constant will be removed after Dune 2.11")]] const static int dofsPerQuad        = (k+1)*(k+1);
  [[deprecated("This constant will be removed after Dune 2.11")]] const static int dofsPerTetrahedron = (k+1)*(k+2)*(k+3)/6;
  [[deprecated("This constant will be removed after Dune 2.11")]] const static int dofsPerPrism       = (k+1)*(k+1)*(k+2)/2;
  [[deprecated("This constant will be removed after Dune 2.11")]] const static int dofsPerHexahedron  = (k+1)*(k+1)*(k+1);
  [[deprecated("This constant will be removed after Dune 2.11")]] const static int dofsPerPyramid     = (k+1)*(k+2)*(2*k+3)/6;

  using Node = LagrangeDGNode<GV, k, R>;

  /**
   * \brief Constructor for a given grid view object
   *
   * This constructor requires that the order is given
   * at compile-time using `k>=0`
   */
  LagrangeDGPreBasis(const GridView& gv)
    requires(k>=0)
    : Base(gv, dofLayout(k))
    , order_(k)
  {}

  /**
   * \brief Constructor for a given grid view object
   *
   * This constructor has to be used to set the order
   * dynamically which is enables using `k<0`.
   */
  LagrangeDGPreBasis(const GridView& gv, unsigned int order)
    requires(k<0)
    : Base(gv, dofLayout(order))
    , order_(order)
  {}

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{order_};
  }

  template<class Node, class It>
  It indices(const Node& node, It it) const
  {
    size_type elementOffset = Base::mapper_.index(node.element());
    for(auto i : Dune::range(node.size()))
    {
      *it = {{ (size_type)(elementOffset+i) }};
      ++it;
    }
    return it;
  }

  //! Polynomial order used in the local Lagrange finite-elements
  unsigned int order() const
  {
    return (useDynamicOrder) ? order_ : k;
  }

protected:

  unsigned int order_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a LagrangeDG pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam order  The polynomial order of the ansatz functions
 * \tparam R      The range type of the local basis
 */
template<std::size_t order, typename R=double>
auto lagrangeDG()
{
  return [](const auto& gridView) {
    return LagrangeDGPreBasis<std::decay_t<decltype(gridView)>, order, R>(gridView);
  };
}

/**
 * \brief Create a pre-basis factory that can create a LagrangeDG pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R      The range type of the local basis
 * \param order   The polynomial order of the ansatz functions
 */
template<typename R=double>
auto lagrangeDG(unsigned int order)
{
  return [order](const auto& gridView) {
    return LagrangeDGPreBasis<std::decay_t<decltype(gridView)>, -1, R>(gridView, order);
  };
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Lagrangean-DG finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 * \tparam R The range type of the local basis
 */
template<typename GV, int k=-1, typename R=double>
using LagrangeDGBasis = DefaultGlobalBasis<LagrangeDGPreBasis<GV, k, R> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
