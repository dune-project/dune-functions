// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>
#include <dune/localfunctions/meta/discontinuous.hh>

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
class LagrangeDGNode :
  public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

  // A simple cache handing storing exactly one LFE.  This can be
  // used for grids with only a single element type. In contrast to
  // StaticLagrangeLocalFiniteElementCache this also supports
  // Lagrange*LocalFiniteElement with run-time order.
  template <class LFE>
  class SingleLocalFiniteElementCache
  {
    LFE lfe_;
  public:
    using FiniteElementType = LFE;

    template<class... Args>
    SingleLocalFiniteElementCache(Args&&... args)
      : lfe_(std::forward<Args>(args)...)
    {}

    //! Obtain the cached local finite-element.
    const FiniteElementType& get ([[maybe_unused]] Dune::GeometryType type) const
    {
      return lfe_;
    }
  };

  template<class LFE>
  using DG = Dune::DiscontinuousLocalFiniteElement<LFE>;

  template<class D, class RR, std::size_t dim, int compileTimeOrder>
  class ImplementedLagrangeDGFiniteElements : public Dune::Impl::ImplementedLagrangeFiniteElements<D, RR, dim, compileTimeOrder>
  {
    using Base = Dune::Impl::ImplementedLagrangeFiniteElements<D, RR, dim, compileTimeOrder>;
  public:

    using Base::Base;

    auto getImplementations() const
    {
      return std::apply([&](auto... indexLFEPair) {
        return std::make_tuple(
          std::make_pair(std::get<0>(indexLFEPair), [lfe = std::get<1>(indexLFEPair)()]() { return DG<std::decay_t<decltype(lfe)>>(lfe); })...
        );
      }, Base::getImplementations());
    }
  };

  template<class D, class RR, std::size_t dim, int compileTimeOrder = -1>
  using LagrangeDGLocalFiniteElementCache = LocalFiniteElementVariantCache<ImplementedLagrangeDGFiniteElements<D,RR,dim, compileTimeOrder>>;

  // Utility function to construct the FiniteElementCache type.
  // Since the function is just a helper to generate a type,
  // it hands out a MetaType<T> instead of a raw T.
  static constexpr auto makeCacheType()
  {
    using D = typename GV::ctype;
    if constexpr (Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>::v)
    {
      constexpr auto type = Dune::GeometryType(Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId, GV::dimension);
      if constexpr (type.isSimplex())
        return Dune::MetaType<SingleLocalFiniteElementCache<DG<Dune::LagrangeSimplexLocalFiniteElement<D,R,dim,k>>>>{};
      else if constexpr (type.isCube())
        return Dune::MetaType<SingleLocalFiniteElementCache<DG<Dune::LagrangeCubeLocalFiniteElement<D,R,dim,k>>>>{};
      else if constexpr (type.isPrism())
        return Dune::MetaType<SingleLocalFiniteElementCache<DG<Dune::LagrangePrismLocalFiniteElement<D,R,k>>>>{};
      else if constexpr (type.isPyramid())
        return Dune::MetaType<SingleLocalFiniteElementCache<DG<Dune::LagrangePyramidLocalFiniteElement<D,R,k>>>>{};
    }
    else
      return Dune::MetaType<LagrangeDGLocalFiniteElementCache<D,R,dim,k>>{};
  }

  using FiniteElementCache = typename decltype(makeCacheType())::type;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  //! Constructor without order (uses the compile-time value)
  LagrangeDGNode() :
    LagrangeDGNode(k)
  {}

  //! Constructor with a run-time order
  LagrangeDGNode(unsigned int order) :
    cache_(order),
    finiteElement_(nullptr),
    element_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
    this->setSize(finiteElement_->size());
  }

protected:

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



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
