// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/raviartthomas.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas0cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas0cube3d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas02d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube3d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas12d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas2cube2d.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>

namespace Dune {
namespace Functions {

namespace Impl {

  template<int dim, typename D, typename R, std::size_t k>
  struct RaviartThomasSimplexLocalInfo
  {
    static_assert((AlwaysFalse<D>::value),"The requested type of Raviart-Thomas element is not implemented, sorry!");
  };

  template<typename D, typename R>
  struct RaviartThomasSimplexLocalInfo<2,D,R,0>
  {
    using FiniteElement = RT02DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 8;
  };

  template<typename D, typename R>
  struct RaviartThomasSimplexLocalInfo<2,D,R,1>
  {
    using FiniteElement = RT12DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 8;
  };

  template<int dim, typename D, typename R, std::size_t k>
  struct RaviartThomasCubeLocalInfo
  {
    static_assert((AlwaysFalse<D>::value),"The requested type of Raviart-Thomas element is not implemented, sorry!");
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<2,D,R,0>
  {
    using FiniteElement = RT0Cube2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 16;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<2,D,R,1>
  {
    using FiniteElement = RT1Cube2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 16;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<2,D,R,2>
  {
    using FiniteElement = RT2Cube2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 16;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<3,D,R,0>
  {
    using FiniteElement = RT0Cube3DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 64;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<3,D,R,1>
  {
    using FiniteElement = RT1Cube3DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 64;
  };

  template<typename GV, int dim, typename R, std::size_t k>
  class RaviartThomasLocalFiniteElementMap
  {
    using D = typename GV::ctype;
    constexpr static bool hasFixedElementType = Capabilities::hasSingleGeometryType<typename GV::Grid>::v;

    using CubeFiniteElement    = typename RaviartThomasCubeLocalInfo<dim, D, R, k>::FiniteElement;
    using SimplexFiniteElement = typename RaviartThomasSimplexLocalInfo<dim, D, R, k>::FiniteElement;
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

    RaviartThomasLocalFiniteElementMap(const GV& gv)
      : is_(&(gv.indexSet())), orient_(gv.size(0))
    {
      cubeVariant_.resize(RaviartThomasCubeLocalInfo<dim, D, R, k>::Variants);
      simplexVariant_.resize(RaviartThomasSimplexLocalInfo<dim, D, R, k>::Variants);

      // create all variants
      for (size_t i = 0; i < cubeVariant_.size(); i++)
        cubeVariant_[i] = std::make_unique<CubeFiniteElementImp>(CubeFiniteElement(i));

      for (size_t i = 0; i < simplexVariant_.size(); i++)
        simplexVariant_[i] = std::make_unique<SimplexFiniteElementImp>(SimplexFiniteElement(i));

      // compute orientation for all elements
      // loop once over the grid
      for(const auto& cell : elements(gv))
      {
        unsigned int myId = is_->index(cell);
        orient_[myId] = 0;

        for (const auto& intersection : intersections(gv,cell))
        {
          if (intersection.neighbor() && (is_->index(intersection.outside()) > myId))
            orient_[myId] |= (1 << intersection.indexInInside());
        }
      }
    }

    /** \brief Get local basis functions for entity if the GeometryType is run-time information
     *
     * The AlwaysTrue class is needed because for SFINAE to work the SFINAE argument has to
     * depend on the template argument.
     */
    template<class EntityType,
             std::enable_if_t<!hasFixedElementType and AlwaysTrue<EntityType>::value,int> = 0>
    const FiniteElement& find(const EntityType& e) const
    {
      if (e.type().isCube())
        return *cubeVariant_[orient_[is_->index(e)]];
      else
        return *simplexVariant_[orient_[is_->index(e)]];
    }

    /** \brief Get local basis functions for entity if the GeometryType is known to be a cube
     *
     * The AlwaysTrue class is needed because for SFINAE to work the SFINAE argument has to
     * depend on the template argument.
     */
    template<class EntityType,
             std::enable_if_t<hasFixedElementType and type.isCube() and AlwaysTrue<EntityType>::value,int> = 0>
    const FiniteElement& find(const EntityType& e) const
    {
      return *cubeVariant_[orient_[is_->index(e)]];
    }

    /** \brief Get local basis functions for entity if the GeometryType is known to be a simplex
     *
     * The AlwaysTrue class is needed because for SFINAE to work the SFINAE argument has to
     * depend on the template argument.
     */
    template<class EntityType,
             std::enable_if_t<hasFixedElementType and type.isSimplex() and AlwaysTrue<EntityType>::value,int> = 0>
    const FiniteElement& find(const EntityType& e) const
    {
      return *simplexVariant_[orient_[is_->index(e)]];
    }

    private:
      std::vector<std::unique_ptr<CubeFiniteElementImp> > cubeVariant_;
      std::vector<std::unique_ptr<SimplexFiniteElementImp> > simplexVariant_;
      const typename GV::IndexSet* is_;
      std::vector<unsigned char> orient_;
  };


} // namespace Impl


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   RaviartThomasPreBasis
//   RaviartThomasNodeIndexSet
//   RaviartThomasNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k>
class RaviartThomasNode;

template<typename GV, int k, class MI>
class RaviartThomasNodeIndexSet;

template<typename GV, int k, class MI>
class RaviartThomasPreBasis
{
  static const int dim = GV::dimension;
  using FiniteElementMap = typename Impl::RaviartThomasLocalFiniteElementMap<GV, dim, double, k>;

  template<typename, int, class>
  friend class RaviartThomasNodeIndexSet;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  using Node = RaviartThomasNode<GV, k>;

  using IndexSet = RaviartThomasNodeIndexSet<GV, k, MI>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  RaviartThomasPreBasis(const GridView& gv) :
    gridView_(gv),
    finiteElementMap_(gv)
  {
    // There is no inherent reason why the basis shouldn't work for grids with more than one
    // element types.  Somebody simply has to sit down and implement the missing bits.
    if (gv.indexSet().types(0).size() > 1)
      DUNE_THROW(Dune::NotImplemented, "Raviart-Thomas basis is only implemented for grids with a single element type");

    GeometryType type = gv.template begin<0>()->type();
    const static int dofsPerElement = (dim == 2) ? (type.isCube() ? k*(k+1)*dim : k*dim) : k*(k+1)*(k+1)*dim;
    constexpr int dofsPerFace    = (dim == 2) ? k+1 : 3*k+1;

    dofsPerCodim_ = {{dofsPerElement, dofsPerFace}};
  }

  void initializeIndices()
  {
    codimOffset_[0] = 0;
    codimOffset_[1] = codimOffset_[0] + dofsPerCodim_[0] * gridView_.size(0);
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
    return dofsPerCodim_[0] * gridView_.size(0) + dofsPerCodim_[1] * gridView_.size(1);
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    size_type result = 0;
    for (auto&& type : gridView_.indexSet().types(0))
    {
      size_t numFaces = ReferenceElements<double,dim>::general(type).size(1);
      const static int dofsPerCodim0 = (dim == 2) ? (type.isCube() ? k*(k+1)*dim : k*dim) : k*(k+1)*(k+1)*dim;
      constexpr int dofsPerCodim1 = (dim == 2) ? k+1 : 3*k+1;
      result = std::max(result, dofsPerCodim0 + dofsPerCodim1 * numFaces);
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



template<typename GV, int k>
class RaviartThomasNode :
  public LeafBasisNode
{
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElementMap = typename Impl::RaviartThomasLocalFiniteElementMap<GV, dim, double, k>;
  using FiniteElement = typename FiniteElementMap::FiniteElement;

  RaviartThomasNode(const FiniteElementMap* finiteElementMap) :
    finiteElement_(nullptr),
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
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(finiteElementMap_->find(*element_));
    this->setSize(finiteElement_->size());
  }

protected:

  const FiniteElement* finiteElement_;
  const Element* element_;
  const FiniteElementMap* finiteElementMap_;
};

template<typename GV, int k, class MI>
class RaviartThomasNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = RaviartThomasPreBasis<GV, k, MI>;

  using Node = RaviartThomasNode<GV,k>;

  RaviartThomasNodeIndexSet(const PreBasis& preBasis) :
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
   *
   * This assumes dim \in \lbrace 2, 3 \rbrace.
   */
  template<typename It>
  It indices(It it) const
  {
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    // throw if Element is not of predefined type
    if (not(element.type().isCube()) and not(element.type().isSimplex()))
      DUNE_THROW(Dune::NotImplemented, "RaviartThomasBasis only implemented for cube and simplex elements.");

    for(std::size_t i=0, end=size(); i<end; ++i, ++it)
    {
      Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);

      // The dimension of the entity that the current dof is related to
      size_t subentity = localKey.subEntity();
      size_t codim = localKey.codim();

      if (not(codim==0 or codim==1))
        DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the RaviartThomasBasis");

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

namespace Imp {

template<std::size_t k>
class RaviartThomasPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return RaviartThomasPreBasis<GridView, k, MultiIndex>(gridView);
  }

};

} // end namespace BasisFactory::Imp

/**
 * \brief Create a pre-basis factory that can create a Raviart-Thomas pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k Order of the Raviart-Thomas element
 */
template<std::size_t k>
auto raviartThomas()
{
  return Imp::RaviartThomasPreBasisFactory<k>();
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a k-th-order Raviart Thomas finite element space
 *
 * TODO: Fix this for grids with more than one element type
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using RaviartThomasBasis = DefaultGlobalBasis<RaviartThomasPreBasis<GV, k, FlatMultiIndex<std::size_t>> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASBASIS_HH
