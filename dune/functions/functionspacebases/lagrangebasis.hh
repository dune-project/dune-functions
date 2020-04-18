// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS2_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS2_HH

#include <type_traits>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

namespace Impl {

  // Cache for LocalFiniteElements for the lagrange pointset
  template <class Domain, class Range, int dim>
  class LagrangeLFECache
  {
  public:
    using FiniteElementType = LagrangeLocalFiniteElement<EquidistantPointSet, dim, Domain, Range>;

    LagrangeLFECache (unsigned int order)
      : order_(order)
    {}

    FiniteElementType const& get (GeometryType type)
    {
      auto it = data_.find(type);
      if (it == data_.end())
        it = data_.emplace(type,FiniteElementType(type,order_)).first;
      return it->second;
    }

  private:
    unsigned int order_;
    std::map<GeometryType, FiniteElementType> data_;
  };

  // Cache for LocalFiniteElements of a fixed GeometryType
  template <GeometryType::Id id, class Domain, class Range, int dim, int k>
  class FixedGeometryTypeLagrangeLFECache
  {
    struct UnknownToplogy {};

    static constexpr bool isSimplex = GeometryType(id).isSimplex();
    static constexpr bool isCube = GeometryType(id).isCube();
    static constexpr bool isPrism = GeometryType(id).isPrism();
    static constexpr bool isPyramid = GeometryType(id).isPyramid();

  public:
    using FiniteElementType
      = std::conditional_t<isSimplex, LagrangeSimplexLocalFiniteElement<Domain,Range,dim,k>,
        std::conditional_t<isCube,    LagrangeCubeLocalFiniteElement<Domain,Range,dim,k>,
        std::conditional_t<isPrism,   LagrangePrismLocalFiniteElement<Domain,Range,k>,
        std::conditional_t<isPyramid, LagrangePyramidLocalFiniteElement<Domain,Range,k>, UnknownToplogy> > > >;

    template <class... Args>
    FiniteElementType const& get (GeometryType type, Args&&... args)
    {
      assert(GeometryType::Id(type) == id);

      // always recreate local finite-element since e.g. vertex map could have changed and
      // different arguments are given
      lfe_.emplace(std::forward<Args>(args)...);
      return *lfe_;
    }

  private:
    std::optional<FiniteElementType> lfe_;
  };


  //! Creator for lagrange LFECache with runtime polynomial order
  template <class Range>
  struct DynamicLagrangeLFECache
  {
    template <class GV>
    using type = LagrangeLFECache<typename GV::ctype,Range,GV::dimension>;
  };

  //! Creator for lagrange LFECache with compiletime polynomial order `k`
  template <class Range, std::size_t k>
  struct StaticLagrangeLFECache
  {
    template <class GV>
    using type = PQkLocalFiniteElementCache<typename GV::ctype, Range, GV::dimension, k>;
  };


  //! Creator for lagrange LFECache for a fixed GeometryType
  template <class Range, std::size_t k>
  struct FixedGeometryTypeLFECache
  {
    template <class GV>
    static constexpr GeometryType::Id id
      = GeometryType(Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId, GV::dimension);

    template <class GV>
    using type = FixedGeometryTypeLagrangeLFECache<id<GV>, typename GV::ctype, Range, GV::dimension, k>;
  };

  //! Select the static LFECache if k > 0, else the dynamic LFECache
  template <class GV, int k, class Range>
  using LFECacheSelector = std::conditional_t<(k >= 0),
    std::conditional_t<Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>::v,
      typename FixedGeometryTypeLFECache<Range,k>::template type<GV>,
      typename StaticLagrangeLFECache<Range,k>::template type<GV>
      >,
    typename DynamicLagrangeLFECache<Range>::template type<GV>
    >;

} // end namespace Impl


// *****************************************************************************
// This is the reusable part of the LagrangeBasis. It contains
//
//   LagrangePreBasis
//   LagrangeNodeIndexSet
//   LagrangeNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template <class GV, class LFECache>
class LagrangeNode;

template <class GV, class LFECache, class MI>
class LagrangeNodeIndexSet;

template <class GV, class LFECache, class MI>
class LagrangePreBasis;



/**
 * \brief A pre-basis for a PQ-lagrange bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV         The grid view that the FE basis is defined on
 * \tparam LFECache   Cache for constructing the LocalFiniteElement type for a GeometryType
 * \tparam MI         Type to be used for multi-indices
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 */
template <class GV, class LFECache, class MI>
class LagrangePreBasis
{
  static const int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

private:

  template <class, class, class>
  friend class LagrangeNodeIndexSet;

public:

  //! Template mapping root tree path to type of created tree node
  using Node = LagrangeNode<GV, LFECache>;

  //! Template mapping root tree path to type of created tree node index set
  using IndexSet = LagrangeNodeIndexSet<GV, LFECache, MI>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object and additional argument to
  //! construct the LFECache.
  template <class... Args,
    std::enable_if_t<std::is_constructible_v<LFECache, Args...>, int> = 0>
  LagrangePreBasis (const GridView& gv, Args&&... args)
    : gridView_(gv)
    , lfeCache_(std::forward<Args>(args)...)
    , order_(lfeCache_.get(GeometryTypes::simplex(dim)).localBasis().order())
  {
    for (int i=0; i<=dim; i++)
    {
      dofsPerCube_[i] = computeDofsPerCube(i);
      dofsPerSimplex_[i] = computeDofsPerSimplex(i);
    }
    dofsPerPrism_ = computeDofsPerPrism();
    dofsPerPyramid_ = computeDofsPerPyramid();
  }

  //! Initialize the global indices
  void initializeIndices ()
  {
    vertexOffset_        = 0;
    edgeOffset_            = vertexOffset_          + dofsPerCube(0) * ((size_type)gridView_.size(dim));

    if (dim>=2)
    {
      triangleOffset_      = edgeOffset_            + dofsPerCube(1) * ((size_type) gridView_.size(dim-1));

      quadrilateralOffset_ = triangleOffset_        + dofsPerSimplex(2) * ((size_type)gridView_.size(Dune::GeometryTypes::triangle));
    }

    if (dim==3) {
      tetrahedronOffset_   = quadrilateralOffset_ + dofsPerCube(2) * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral));

      prismOffset_         = tetrahedronOffset_   +   dofsPerSimplex(3) * ((size_type)gridView_.size(Dune::GeometryTypes::tetrahedron));

      hexahedronOffset_    = prismOffset_         +   dofsPerPrism() * ((size_type)gridView_.size(Dune::GeometryTypes::prism));

      pyramidOffset_       = hexahedronOffset_    +   dofsPerCube(3) * ((size_type)gridView_.size(Dune::GeometryTypes::hexahedron));
    }
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView () const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  /**
   * \brief Create tree node
   */
  Node makeNode () const
  {
    return Node{lfeCache_};
  }

  /**
   * \brief Create tree node index set
   *
   * Create an index set suitable for the tree node obtained
   * by makeNode().
   */
  IndexSet makeIndexSet () const
  {
    return IndexSet{*this};
  }

  //! Same as size(prefix) with empty prefix
  size_type size () const
  {
    switch (dim)
    {
      case 1:
        return dofsPerCube(0) * ((size_type)gridView_.size(dim))
          + dofsPerCube(1) * ((size_type)gridView_.size(dim-1));
      case 2:
      {
        return dofsPerCube(0) * ((size_type)gridView_.size(dim))
          + dofsPerCube(1) * ((size_type)gridView_.size(dim-1))
          + dofsPerSimplex(2) * ((size_type)gridView_.size(Dune::GeometryTypes::triangle))
          + dofsPerCube(2) * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral));
      }
      case 3:
      {
        return dofsPerCube(0) * ((size_type)gridView_.size(dim))
          + dofsPerCube(1) * ((size_type)gridView_.size(dim-1))
          + dofsPerSimplex(2) * ((size_type)gridView_.size(Dune::GeometryTypes::triangle))
          + dofsPerCube(2) * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral))
          + dofsPerSimplex(3) * ((size_type)gridView_.size(Dune::GeometryTypes::tetrahedron))
          + dofsPerPyramid() * ((size_type)gridView_.size(Dune::GeometryTypes::pyramid))
          + dofsPerPrism() * ((size_type)gridView_.size(Dune::GeometryTypes::prism))
          + dofsPerCube(3) * ((size_type)gridView_.size(Dune::GeometryTypes::hexahedron));
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Return number of possible values for next position in multi index
  size_type size (const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension () const
  {
    return size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize () const
  {
    // That cast to unsigned int is necessary because GV::dimension is an enum,
    // which is not recognized by the power method as an integer type...
    return power(order()+1, (unsigned int)GV::dimension);
  }

protected:
  unsigned int order () const
  {
    return order_ ;
  }

  //! Number of degrees of freedom assigned to a simplex (without the ones assigned to its faces!)
  size_type dofsPerSimplex (std::size_t simplexDim) const
  {
    return dofsPerSimplex_[simplexDim];
  }

  //! Number of degrees of freedom assigned to a cube (without the ones assigned to its faces!)
  size_type dofsPerCube (std::size_t cubeDim) const
  {
    return dofsPerCube_[cubeDim];
  }

  size_type dofsPerPrism () const
  {
    return dofsPerPrism_;
  }

  size_type dofsPerPyramid () const
  {
    return dofsPerPyramid_;
  }

  //! Number of degrees of freedom assigned to a simplex (without the ones assigned to its faces!)
  size_type computeDofsPerSimplex (std::size_t simplexDim) const
  {
    return order() == 0 ? (dim == simplexDim ? 1 : 0) : Dune::binomial(std::size_t(order()-1),simplexDim);
  }

  //! Number of degrees of freedom assigned to a cube (without the ones assigned to its faces!)
  size_type computeDofsPerCube (std::size_t cubeDim) const
  {
    return order() == 0 ? (dim == cubeDim ? 1 : 0) : Dune::power(order()-1, cubeDim);
  }

  size_type computeDofsPerPrism () const
  {
    return order() == 0 ? (dim == 3 ? 1 : 0) : (order()-1)*(order()-1)*(order()-2)/2;
  }

  size_type computeDofsPerPyramid () const
  {
    return order() == 0 ? (dim == 3 ? 1 : 0) : (order()-2)*(order()-1)*(2*order()-3)/6;
  }

protected:
  GridView gridView_;
  LFECache lfeCache_;
  unsigned int order_;

  std::array<size_type,dim+1> dofsPerSimplex_;
  std::array<size_type,dim+1> dofsPerCube_;
  size_type dofsPerPrism_;
  size_type dofsPerPyramid_;

  size_type vertexOffset_;
  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type quadrilateralOffset_;
  size_type tetrahedronOffset_;
  size_type pyramidOffset_;
  size_type prismOffset_;
  size_type hexahedronOffset_;
};

template <class GV, int k, class MI, class R = double>
using StaticLagrangePreBasis = LagrangePreBasis<GV,Impl::LFECacheSelector<GV,k,R>,MI>;


template <class GV, class LFECache>
class LagrangeNode
    : public LeafBasisNode
{
public:
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename LFECache::FiniteElementType;

public:
  //! Constructor that creates the local-FE cache
  template <class... Args,
    std::enable_if_t<std::is_constructible_v<LFECache, Args...>, int> = 0>
  LagrangeNode (Args&&... args)
    : lfeCache_(std::forward<Args>(args)...)
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
    finiteElement_ = &(lfeCache_.get(element_->type()));
    this->setSize(finiteElement_->size());
  }

protected:
  LFECache lfeCache_;

  const FiniteElement* finiteElement_ = nullptr;
  const Element* element_ = nullptr;
};

template <class GV, int k, class R = double>
using StaticLagrangeNode = LagrangeNode<GV,Impl::LFECacheSelector<GV,k,R>>;


template <class GV, class LFECache, class MI>
class LagrangeNodeIndexSet
{
  enum { dim = GV::dimension };

public:
  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = LagrangePreBasis<GV, LFECache, MI>;
  using Node = LagrangeNode<GV, LFECache>;

  LagrangeNodeIndexSet (const PreBasis& preBasis)
    : preBasis_(&preBasis)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind (const Node& node)
  {
    node_ = &node;
  }

  /** \brief Unbind the view
   */
  void unbind ()
  {
    node_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size () const
  {
    assert(node_ != nullptr);
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template <class It>
  It indices (It it) const
  {
    assert(node_ != nullptr);
    for (size_type i = 0, end = node_->finiteElement().size() ; i < end ; ++it, ++i)
      {
        Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
        const auto& gridIndexSet = preBasis_->gridView().indexSet();
        const auto& element = node_->element();

        // The dimension of the entity that the current dof is related to
        auto dofDim = dim - localKey.codim();

        // Test for a vertex dof
        // The test for k==1 is redundant, but having it here allows the compiler to conclude
        // at compile-time that the dofDim==0 case is the only one that will ever happen.
        // This leads to measurable speed-up: see
        //   https://gitlab.dune-project.org/staging/dune-functions/issues/30
        if (order()==1 || dofDim==0) {
          *it = {{ (size_type)(gridIndexSet.subIndex(element,localKey.subEntity(),dim)) }};
          continue;
        }

        if (dofDim==1)
          {  // edge dof
            if (dim==1)  // element dof -- any local numbering is fine
              {
                *it = {{ preBasis_->edgeOffset_
                         + preBasis_->dofsPerCube(1) * ((size_type)gridIndexSet.subIndex(element,0,0))
                         + localKey.index() }};
                continue;
              }
            else
              {
                const auto refElement
                  = Dune::referenceElement<double,dim>(element.type());

                // we have to reverse the numbering if the local triangle edge is
                // not aligned with the global edge
                auto v0 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
                auto v1 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
                bool flip = (v0 > v1);
                *it = {{ (flip)
                         ? preBasis_->edgeOffset_
                         + preBasis_->dofsPerCube(1)*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                         + (preBasis_->dofsPerCube(1)-1)-localKey.index()
                         : preBasis_->edgeOffset_
                         + preBasis_->dofsPerCube(1)*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                         + localKey.index() }};
                continue;
              }
          }

        if (dofDim==2)
          {
            if (dim==2)   // element dof -- any local numbering is fine
              {
                if (element.type().isTriangle())
                  {
                    *it = {{ preBasis_->triangleOffset_ + preBasis_->dofsPerSimplex(2)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isQuadrilateral())
                  {
                    *it = {{ preBasis_->quadrilateralOffset_ + preBasis_->dofsPerCube(2)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else
                  DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
              } else
              {
                const auto refElement
                  = Dune::referenceElement<double,dim>(element.type());

                if (order()>3)
                  DUNE_THROW(Dune::NotImplemented, "LagrangeNodalBasis for 3D grids is only implemented if k<=3");

                if (order()==3 and !refElement.type(localKey.subEntity(), localKey.codim()).isTriangle())
                  DUNE_THROW(Dune::NotImplemented, "LagrangeNodalBasis for 3D grids with k==3 is only implemented if the grid is a simplex grid");

                *it = {{ preBasis_->triangleOffset_ + ((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())) }};
                continue;
              }
          }

        if (dofDim==3)
          {
            if (dim==3)   // element dof -- any local numbering is fine
              {
                if (element.type().isTetrahedron())
                  {
                    *it = {{ preBasis_->tetrahedronOffset_ + preBasis_->dofsPerSimplex(3)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isHexahedron())
                  {
                    *it = {{ preBasis_->hexahedronOffset_ + preBasis_->dofsPerCube(3)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isPrism())
                  {
                    *it = {{ preBasis_->prismOffset_ + preBasis_->dofsPerPrism()*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isPyramid())
                  {
                    *it = {{ preBasis_->pyramidOffset_ + preBasis_->dofsPerPyramid()*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else
                  DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedra, hexahedra, prisms, or pyramids");
              } else
              DUNE_THROW(Dune::NotImplemented, "Grids of dimension larger than 3 are no supported");
          }
        DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the LagrangeNodalBasis");
      }
    return it;
  }

protected:
  unsigned int order() const
  {
    return preBasis_->order();
  }

protected:
  const PreBasis* preBasis_;
  const Node* node_ = nullptr;
};


namespace BasisFactory {
namespace Impl {

template <template <class> class LFECache, bool useDynamicOrder = false>
class LagrangePreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  // \brief Constructor for factory with compile-time order
  LagrangePreBasisFactory ()
    : order_(0)
  {}

  // \brief Constructor for factory with run-time order (template argument k is disregarded)
  LagrangePreBasisFactory (unsigned int order)
    : order_(order)
  {}

  template <class MultiIndex, class GridView>
  auto makePreBasis (const GridView& gridView) const
  {
    using Cache = LFECache<GridView>;

    if constexpr (useDynamicOrder)
      return LagrangePreBasis<GridView, Cache, MultiIndex>(gridView, order_);
    else
      return LagrangePreBasis<GridView, Cache, MultiIndex>(gridView);
  }

private:
  unsigned int order_;
};

} // end namespace BasisFactory::Impl


/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions
 * \tparam R   The range type of the local basis
 */
template <std::size_t k, class R = double>
auto lagrange ()
{
  return Impl::LagrangePreBasisFactory<Functions::Impl::StaticLagrangeLFECache<R,k>::template type>();
}

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis with a run-time order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 */
template <class R = double>
auto lagrange (int order)
{
  return Impl::LagrangePreBasisFactory<Functions::Impl::DynamicLagrangeLFECache<R>::template type,true>(order);
}

} // end namespace BasisFactory


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
 * of LagrangePreBasis.
 *
 * \warning The implementation of the basis with run-time order order uses the
 *   LagrangeFiniteElement implementation of dune-localfunctions, which is known
 *   to violate strict-aliasing rules
 *   (see https://gitlab.dune-project.org/core/dune-localfunctions/issues/14)
 *   Keep this in mind if ever you experience difficult-to-explain crashes
 *   or wrong results.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis; -1 means 'order determined at run-time'
 * \tparam R The range type of the local basis
 */
template <class GV, int k = -1, class R = double>
using LagrangeBasis
  = DefaultGlobalBasis<LagrangePreBasis<GV, Impl::LFECacheSelector<GV,k,R>, FlatMultiIndex<std::size_t>> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS2_HH
