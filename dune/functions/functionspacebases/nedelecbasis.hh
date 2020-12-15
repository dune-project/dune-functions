// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NEDELECBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NEDELECBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>
#include <dune/localfunctions/nedelec.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
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

    // TODO: The cube code is an investment for the future.  We don't have local finite elements for cubes yet,
    //       but let's have the infrastructure for them already.
    using CubeFiniteElement    = Nedelec1stKindSimplexLocalFiniteElement<D,R,dim,order>;
    using SimplexFiniteElement = Nedelec1stKindSimplexLocalFiniteElement<D,R,dim,order>;
    using CubeFiniteElementImp = typename std::conditional<hasFixedElementType,
                                                           CubeFiniteElement,
                                                           LocalFiniteElementVirtualImp<CubeFiniteElement> >::type;
    using SimplexFiniteElementImp = typename std::conditional<hasFixedElementType,
                                                              SimplexFiniteElement,
                                                              LocalFiniteElementVirtualImp<SimplexFiniteElement> >::type;

  public:

    using T = LocalBasisTraits<D, dim, FieldVector<D,dim>, R, dim, FieldVector<R,dim>, FieldMatrix<D,dim,dim> >;

    constexpr static unsigned int  topologyId = Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId;  // meaningless if hasFixedElementType is false
    constexpr static GeometryType type = GeometryType(topologyId, GV::dimension);

    using FiniteElement = typename std::conditional<hasFixedElementType,
                                           typename std::conditional<type.isCube(),CubeFiniteElement,SimplexFiniteElement>::type,
                                           LocalFiniteElementVirtualInterface<T> >::type;

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
      cubeVariant_.resize(numVariants(GeometryTypes::cube(dim)));
      simplexVariant_.resize(numVariants(GeometryTypes::simplex(dim)));

      // create all variants
      for (size_t i = 0; i < cubeVariant_.size(); i++)
        cubeVariant_[i] = std::make_unique<CubeFiniteElementImp>(CubeFiniteElement(i));

      for (size_t i = 0; i < simplexVariant_.size(); i++)
        simplexVariant_[i] = std::make_unique<SimplexFiniteElementImp>(SimplexFiniteElement(i));

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
      }
    }

    template<class Element>
    const auto& find(const Element& element) const
    {
      if constexpr (!hasFixedElementType)
      {
        return (element.type().isCube()) ? *cubeVariant_[orientation_[elementMapper_.index(element)]]
                                         : *simplexVariant_[orientation_[elementMapper_.index(element)]];
      }
      else
      {
        if constexpr (type.isCube())
          return *cubeVariant_[orientation_[elementMapper_.index(element)]];
        else
          return *simplexVariant_[orientation_[elementMapper_.index(element)]];
      }
    }

    private:
      std::vector<std::unique_ptr<CubeFiniteElementImp> > cubeVariant_;
      std::vector<std::unique_ptr<SimplexFiniteElementImp> > simplexVariant_;
      const Dune::MultipleCodimMultipleGeomTypeMapper<GV> elementMapper_;
      std::vector<unsigned char> orientation_;
  };


} // namespace Impl


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   NedelecPreBasis
//   NedelecNodeIndexSet
//   NedelecNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, typename Range, std::size_t kind, int order>
class NedelecNode;

template<typename GV, typename Range, std::size_t kind, int order, class MI>
class NedelecNodeIndexSet;

template<typename GV, typename Range, std::size_t kind, int order, class MI>
class NedelecPreBasis
{
  static const int dim = GV::dimension;
  static_assert(kind==1, "Only the Nedelec basis of the first kind is currently implemented!");
  using FiniteElementMap = typename Impl::Nedelec1stKindLocalFiniteElementMap<GV, dim, Range, order>;

  friend class NedelecNodeIndexSet<GV,Range,kind,order,MI>;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  using Node = NedelecNode<GV, Range, kind, order>;

  using IndexSet = NedelecNodeIndexSet<GV, Range, kind, order, MI>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  NedelecPreBasis(const GridView& gv) :
    gridView_(gv),
    finiteElementMap_(gv)
  {
    if (kind!=1)
      DUNE_THROW(NotImplemented, "Only Nedelec elements of the first kind are implemented!");

    // There is no inherent reason why the basis shouldn't work for grids with more than one
    // element types.  Somebody simply has to sit down and implement the missing bits.
    if (gv.indexSet().types(0).size() > 1)
      DUNE_THROW(NotImplemented, "Nedelec basis is only implemented for grids with a single element type");

    if (!gv.indexSet().types(0)[0].isSimplex())
      DUNE_THROW(NotImplemented, "Nédélec basis is only implemented for grids with simplex elements.");

    if (order>1)
      DUNE_THROW(NotImplemented, "Only first-order elements are implemented");

    if (dim!=2 && dim!=3)
      DUNE_THROW(NotImplemented, "Only 2d and 3d Nédélec elements are implemented");

    std::fill(dofsPerCodim_.begin(), dofsPerCodim_.end(), 0);

    if (kind==1)
      dofsPerCodim_[dim-1] = 1;  // First-order: One dof per edge, and no others
  }

  void initializeIndices()
  {
    codimOffset_[0] = 0;

    for (std::size_t i=0; i<dim; i++)
      codimOffset_[i+1] = codimOffset_[i] + dofsPerCodim_[i] * gridView_.size(i);
  }

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
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{&finiteElementMap_};
  }

  /**
   * \brief Create tree node index set
   *
   * Create an index set suitable for the tree node obtained
   * by makeNode().
   */
  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }

  size_type size() const
  {
    size_type result = 0;
    for (int i=0; i<=dim; i++)
      result += dofsPerCodim_[i] * gridView_.size(i);
    return result;
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  size_type dimension() const
  {
    return size();
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

protected:
  const GridView gridView_;
  std::array<size_t,dim+1> codimOffset_;
  FiniteElementMap finiteElementMap_;
  // Number of dofs per entity type depending on the entity's codimension and type
  std::array<int,dim+1> dofsPerCodim_;
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

template<typename GV, typename Range, std::size_t kind, int order, class MI>
class NedelecNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = NedelecPreBasis<GV, Range, kind, order, MI>;

  using Node = NedelecNode<GV, Range, kind, order>;

  NedelecNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    node_ = &node;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return node_->finiteElement().size();
  }

  /**
   * \brief Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
   */
  template<typename It>
  It indices(It it) const
  {
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    // throw if Element is not of predefined type
    if (not(element.type().isCube()) and not(element.type().isSimplex()))
      DUNE_THROW(NotImplemented, "NedelecBasis only implemented for cube and simplex elements.");

    for(std::size_t i=0, end=size(); i<end; ++i, ++it)
    {
      Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);

      // The dimension of the entity that the current dof is related to
      size_t subentity = localKey.subEntity();
      size_t codim = localKey.codim();

      *it = { preBasis_->codimOffset_[codim] +
        preBasis_->dofsPerCodim_[codim] * gridIndexSet.subIndex(element, subentity, codim) + localKey.index() };
    }

    return it;
  }

protected:
  const PreBasis* preBasis_;
  const Node* node_;
};



namespace BasisFactory {

namespace Impl {

template<std::size_t kind, std::size_t order, typename Range>
class NedelecPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return NedelecPreBasis<GridView, Range, kind, order, MultiIndex>(gridView);
  }

};

} // end namespace BasisFactory::Impl

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
  return Impl::NedelecPreBasisFactory<kind, order, Range>();
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
using NedelecBasis = DefaultGlobalBasis<NedelecPreBasis<GV, Range, kind, order, FlatMultiIndex<std::size_t> > >;

} // end namespace Dune::Functions


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NEDELECBASIS_HH
