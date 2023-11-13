// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NEDELECBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NEDELECBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localfiniteelementvariant.hh>
#include <dune/localfunctions/nedelec.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
#include <dune/functions/functionspacebases/leafprebasismixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>

namespace Dune::Functions
{

namespace Impl
{
  template<typename GV, int dim, typename R, std::size_t order>
  class Nedelec1stKindLocalFiniteElementMap
  {
    using D = typename GV::ctype;
    constexpr static bool hasFixedElementType = Capabilities::hasSingleGeometryType<typename GV::Grid>::v;

    using CubeFiniteElement    = Nedelec1stKindCubeLocalFiniteElement<D,R,dim,order>;
    using SimplexFiniteElement = Nedelec1stKindSimplexLocalFiniteElement<D,R,dim,order>;

  public:

    using T = LocalBasisTraits<D, dim, FieldVector<D,dim>, R, dim, FieldVector<R,dim>, FieldMatrix<D,dim,dim> >;

    constexpr static unsigned int  topologyId = Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId;  // meaningless if hasFixedElementType is false
    constexpr static GeometryType type = GeometryType(topologyId, GV::dimension);

    using FiniteElement = std::conditional_t<hasFixedElementType,
                                           std::conditional_t<type.isCube(),CubeFiniteElement,SimplexFiniteElement>,
                                           LocalFiniteElementVariant<CubeFiniteElement, SimplexFiniteElement> >;

    static std::size_t numVariants(GeometryType type)
    {
      if (order!=1)  // I am not sure whether the formula below is correct for all orders.
        DUNE_THROW(NotImplemented, "Only Nedelec elements of order 1 are implemented!");

      auto numEdges = referenceElement<D,dim>(type).size(dim-1);
      return power(2,numEdges);
    }

    Nedelec1stKindLocalFiniteElementMap(const GV& gv)
      : elementMapper_(gv, mcmgElementLayout()),
        orientation_(gv.size(0))
    {
      // create all variants
      if constexpr (hasFixedElementType)
      {
        variants_.resize(numVariants(type));
        for (size_t i = 0; i < numVariants(type); i++)
          variants_[i] = FiniteElement(i);
      }
      else
      {
        // for mixed grids add offset for cubes
        variants_.resize(numVariants(GeometryTypes::simplex(dim)) + numVariants(GeometryTypes::cube(dim)));
        for (size_t i = 0; i < numVariants(GeometryTypes::simplex(dim)); i++)
          variants_[i] = SimplexFiniteElement(i);
        for (size_t i = 0; i < numVariants(GeometryTypes::cube(dim)); i++)
          variants_[i + numVariants(GeometryTypes::simplex(dim))] = CubeFiniteElement(i);
      }


      // compute orientation for all elements
      const auto& indexSet = gv.indexSet();

      for(const auto& element : elements(gv))
      {
        const auto& refElement = referenceElement(element);
        auto elementIndex = elementMapper_.index(element);
        orientation_[elementIndex] = 0;

        for (std::size_t i=0; i<element.subEntities(dim-1); i++)
        {
          // Local vertex indices within the element
          auto localV0 = refElement.subEntity(i,dim-1, 0,dim);
          auto localV1 = refElement.subEntity(i,dim-1, 1,dim);

          // Global vertex indices within the grid
          auto globalV0 = indexSet.subIndex(element,localV0,dim);
          auto globalV1 = indexSet.subIndex(element,localV1,dim);

          if ( (localV0<localV1 && globalV0>globalV1) || (localV0>localV1 && globalV0<globalV1) )
            orientation_[elementIndex] |= (1 << i);
        }
        // for mixed grids add offset for cubes
        if constexpr (!hasFixedElementType)
          if (element.type().isCube())
            orientation_[elementIndex] += numVariants(GeometryTypes::simplex(dim));
      }
    }

    template<class Element>
    const auto& find(const Element& element) const
    {
      return variants_[orientation_[elementMapper_.index(element)]];
    }

    private:
      std::vector<FiniteElement> variants_;
      Dune::MultipleCodimMultipleGeomTypeMapper<GV> elementMapper_;
      std::vector<unsigned short> orientation_;
  };


} // namespace Impl


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   NedelecPreBasis
//   NedelecNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, typename Range, std::size_t kind, int order>
class NedelecNode;

template<typename GV, typename Range, std::size_t kind, int order>
class NedelecPreBasis :
  public LeafPreBasisMixin< NedelecPreBasis<GV,Range,kind,order> >
{
  static const int dim = GV::dimension;
  static_assert(kind==1, "Only the Nedelec basis of the first kind is currently implemented!");
  using FiniteElementMap = typename Impl::Nedelec1stKindLocalFiniteElementMap<GV, dim, Range, order>;

  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  using Node = NedelecNode<GV, Range, kind, order>;

  /** \brief Constructor for a given grid view object */
  NedelecPreBasis(const GridView& gv) :
    gridView_(gv),
    finiteElementMap_(gv),
    mapper_(gridView_, mcmgLayout(Dim<1>{}))
  {
    if (kind!=1)
      DUNE_THROW(NotImplemented, "Only Nedelec elements of the first kind are implemented!");

    // There is no inherent reason why the basis shouldn't work for grids with more than two
    // element types.  Somebody simply has to sit down and implement the missing bits.
    if (gv.indexSet().types(0).size() > 2)
      DUNE_THROW(NotImplemented, "Nédélec basis is only implemented for grids with simplex and cube elements");

    for(auto type : gv.indexSet().types(0))
      if (!type.isSimplex() && !type.isCube())
        DUNE_THROW(NotImplemented, "Nédélec basis is only implemented for grids with simplex or cube elements.");

    if (order>1)
      DUNE_THROW(NotImplemented, "Only first-order elements are implemented");

    if (dim!=2 && dim!=3)
      DUNE_THROW(NotImplemented, "Only 2d and 3d Nédélec elements are implemented");
  }

  void initializeIndices()
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  /* \brief Update the stored grid view, to be called if the grid has changed */
  void update (const GridView& gv)
  {
    gridView_ = gv;
    mapper_.update(gridView_);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{&finiteElementMap_};
  }

  size_type dimension() const
  {
    return mapper_.size();
  }

  size_type maxNodeSize() const
  {
    size_type result = 0;
    for (auto&& type : gridView_.indexSet().types(0))
    {
      size_type numEdges = referenceElement<typename GV::ctype,dim>(type).size(dim-1);
      result = std::max(result, numEdges);
    }

    return result;
  }

  /**
   * \brief Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
   */
  template<typename It>
  It indices(const Node& node, It it) const
  {
    const auto& element = node.element();

    // throw if Element is not of predefined type
    if (not(element.type().isCube()) and not(element.type().isSimplex()))
      DUNE_THROW(NotImplemented, "NedelecBasis only implemented for cube and simplex elements.");

    for(std::size_t i=0, end=node.size(); i<end; ++i, ++it)
    {
      Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
      *it = { mapper_.subIndex(element, localKey.subEntity(), localKey.codim()) + localKey.index() };
    }

    return it;
  }

protected:
  GridView gridView_;
  FiniteElementMap finiteElementMap_;
  Mapper mapper_;
};



template<typename GV, typename Range, size_t kind, int order>
class NedelecNode :
  public LeafBasisNode
{
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  static_assert(kind==1, "Only Nedelec elements of the first kind are implemented!");
  using FiniteElementMap = typename Impl::Nedelec1stKindLocalFiniteElementMap<GV, dim, Range, order>;
  using FiniteElement = Impl::GlobalValuedLocalFiniteElement<Impl::CovariantPiolaTransformator,
                                                             typename FiniteElementMap::FiniteElement,
                                                             Element>;

  NedelecNode(const FiniteElementMap* finiteElementMap) :
    element_(nullptr),
    finiteElementMap_(finiteElementMap)
  { }

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
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_.bind((finiteElementMap_->find(*element_)), e);
    this->setSize(finiteElement_.size());
  }

protected:

  FiniteElement finiteElement_;
  const Element* element_;
  const FiniteElementMap* finiteElementMap_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a Nédélec pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam kind Kind of the Nédélec element (1 or 2)
 * \tparam order Order of the Nédélec element (lowest order is '1')
 * \tparam Range Number type used for shape function values
 */
template<std::size_t kind, std::size_t order, typename Range=double>
auto nedelec()
{
  return [](const auto& gridView) {
    return NedelecPreBasis<std::decay_t<decltype(gridView)>, Range, kind, order>(gridView);
  };
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a k-th-order Nédélec finite element space
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam Range Number type used for shape function values
 * \tparam kind Kind of the basis: 1 (for Nédélec element of the first kind) or 2
 * \tparam order The order of the basis (lowest order is '1')
 */
template<typename GV, std::size_t kind, std::size_t order, typename Range=double>
using NedelecBasis = DefaultGlobalBasis<NedelecPreBasis<GV, Range, kind, order > >;

} // end namespace Dune::Functions


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NEDELECBASIS_HH
