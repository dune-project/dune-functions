// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

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

template<typename GV, int k>
class LagrangeNode;

template<typename GV, int k, class MI>
class LagrangeNodeIndexSet;

template<typename GV, int k, class MI>
class LagrangePreBasis;



/**
 * \brief A pre-basis for a PQ-lagrange bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions; -1 means 'order determined at run-time'
 * \tparam MI  Type to be used for multi-indices
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 */
template<typename GV, int k, class MI>
class LagrangePreBasis
{
  static const int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

private:

  template<typename, int, class>
  friend class LagrangeNodeIndexSet;

public:

  //! Template mapping root tree path to type of created tree node
  using Node = LagrangeNode<GV, k>;

  //! Template mapping root tree path to type of created tree node index set
  using IndexSet = LagrangeNodeIndexSet<GV, k, MI>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object with compile-time order
  LagrangePreBasis(const GridView& gv) :
    gridView_(gv)
  {}

  //! Constructor for a given grid view object and run-time order
  LagrangePreBasis(const GridView& gv, unsigned int order) :
    gridView_(gv), order_(order)
  {
    if (k>=0)
      DUNE_THROW(RangeError, "Template argument k has to be -1 when supplying a run-time order!");
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_          + dofsPerVertex * ((size_type)gridView_.size(dim));

    if (dim>=2)
    {
      triangleOffset_      = edgeOffset_            + dofsPerEdge * ((size_type) gridView_.size(dim-1));

      quadrilateralOffset_ = triangleOffset_        + dofsPerTriangle * ((size_type)gridView_.size(Dune::GeometryTypes::triangle));
    }

    if (dim==3) {
      tetrahedronOffset_   = quadrilateralOffset_ + dofsPerQuad * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral));

      prismOffset_         = tetrahedronOffset_   +   dofsPerTetrahedron * ((size_type)gridView_.size(Dune::GeometryTypes::tetrahedron));

      hexahedronOffset_    = prismOffset_         +   dofsPerPrism * ((size_type)gridView_.size(Dune::GeometryTypes::prism));

      pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * ((size_type)gridView_.size(Dune::GeometryTypes::hexahedron));
    }
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
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
  Node makeNode() const
  {
    return Node{order_};
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

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    switch (dim)
    {
      case 1:
        return dofsPerVertex * ((size_type)gridView_.size(dim))
          + dofsPerEdge*((size_type)gridView_.size(dim-1));
      case 2:
      {
        return dofsPerVertex * ((size_type)gridView_.size(dim))
          + dofsPerEdge * ((size_type)gridView_.size(dim-1))
          + dofsPerTriangle * ((size_type)gridView_.size(Dune::GeometryTypes::triangle))
          + dofsPerQuad * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral));
      }
      case 3:
      {
        return dofsPerVertex * ((size_type)gridView_.size(dim))
          + dofsPerEdge * ((size_type)gridView_.size(dim-1))
          + dofsPerTriangle * ((size_type)gridView_.size(Dune::GeometryTypes::triangle))
          + dofsPerQuad * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral))
          + dofsPerTetrahedron * ((size_type)gridView_.size(Dune::GeometryTypes::tetrahedron))
          + dofsPerPyramid * ((size_type)gridView_.size(Dune::GeometryTypes::pyramid))
          + dofsPerPrism * ((size_type)gridView_.size(Dune::GeometryTypes::prism))
          + dofsPerHexahedron * ((size_type)gridView_.size(Dune::GeometryTypes::hexahedron));
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return (k>=0)
      ? StaticPower<(k+1),GV::dimension>::power
      : Power<GV::dimension>::eval(order()+1);
  }

protected:
  GridView gridView_;

  unsigned int order() const
  {
    return (k>=0) ? k : order_;
  }

  // Run-time order, only valid if k<0
  unsigned int order_;

    // Compute the number of dofs per entity type
  const size_type dofsPerVertex =
      order() == 0 ? (dim == 0 ? 1 : 0) : 1;
  const size_type dofsPerEdge =
      order() == 0 ? (dim == 1 ? 1 : 0) : order()-1;
  const size_type dofsPerTriangle =
      order() == 0 ? (dim == 2 ? 1 : 0) : (order()-1)*(order()-2)/2;
  const size_type dofsPerQuad =
      order() == 0 ? (dim == 2 ? 1 : 0) : (order()-1)*(order()-1);
  const size_type dofsPerTetrahedron =
      order() == 0 ? (dim == 3 ? 1 : 0) : (order()-3)*(order()-2)*(order()-1)/6;
  const size_type dofsPerPrism =
      order() == 0 ? (dim == 3 ? 1 : 0) : (order()-1)*(order()-1)*(order()-2)/2;
  const size_type dofsPerHexahedron =
      order() == 0 ? (dim == 3 ? 1 : 0) : (order()-1)*(order()-1)*(order()-1);
  const size_type dofsPerPyramid =
      order() == 0 ? (dim == 3 ? 1 : 0) : (order()-2)*(order()-1)*(2*order()-3)/6;

  size_type vertexOffset_;
  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type quadrilateralOffset_;
  size_type tetrahedronOffset_;
  size_type pyramidOffset_;
  size_type prismOffset_;
  size_type hexahedronOffset_;

};



template<typename GV, int k>
class LagrangeNode :
  public LeafBasisNode
{
  // Stores LocalFiniteElement implementations with run-time order as a function of GeometryType
  template<typename D, typename R, int dim>
  class LagrangeRunTimeLFECache
  {
  public:
    using FiniteElementType = LagrangeLocalFiniteElement<EquidistantPointSet,dim,D,R>;

    const FiniteElementType& get(GeometryType type)
    {
      auto i = data_.find(type);
      if (i==data_.end())
        i = data_.emplace(type,FiniteElementType(type,order_)).first;
      return (*i).second;
    }

    std::map<GeometryType, FiniteElementType> data_;
    unsigned int order_;
  };

  static const int dim = GV::dimension;

  using FiniteElementCache = typename std::conditional<(k>=0),
                                                       PQkLocalFiniteElementCache<typename GV::ctype, double, dim, k>,
                                                       LagrangeRunTimeLFECache<typename GV::ctype, double, dim>
                                                      >::type;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  LagrangeNode(unsigned int order) :
    order_(order),
    finiteElement_(nullptr),
    element_(nullptr)
  {
    // Only the cache for the run-time-order case (i.e., k<0), has the 'order_' member
    Hybrid::ifElse(Std::bool_constant<(k<0)>(),
      [&](auto id) {
        id(cache_).order_ = order;
      },
      [&](auto id) {}
      );
  }

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

  unsigned int order() const
  {
    return (k>=0) ? k : order_;
  }

  // Run-time order, only valid if k<0
  unsigned int order_;

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



template<typename GV, int k, class MI>
class LagrangeNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = LagrangePreBasis<GV, k, MI>;

  using Node = LagrangeNode<GV, k>;

  LagrangeNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis),
    node_(nullptr)
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
    assert(node_ != nullptr);
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(It it) const
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
        if (k==1 || dofDim==0) {
          *it = {{ (size_type)(gridIndexSet.subIndex(element,localKey.subEntity(),dim)) }};
          continue;
        }

        if (dofDim==1)
          {  // edge dof
            if (dim==1)  // element dof -- any local numbering is fine
              {
                *it = {{ preBasis_->edgeOffset_
                         + preBasis_->dofsPerEdge * ((size_type)gridIndexSet.subIndex(element,0,0))
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
                         + preBasis_->dofsPerEdge*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                         + (preBasis_->dofsPerEdge-1)-localKey.index()
                         : preBasis_->edgeOffset_
                         + preBasis_->dofsPerEdge*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
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
                    const int interiorLagrangeNodesPerTriangle = (order()-1)*(order()-2)/2;
                    *it = {{ preBasis_->triangleOffset_ + interiorLagrangeNodesPerTriangle*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isQuadrilateral())
                  {
                    const int interiorLagrangeNodesPerQuadrilateral = (order()-1)*(order()-1);
                    *it = {{ preBasis_->quadrilateralOffset_ + interiorLagrangeNodesPerQuadrilateral*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
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
                    *it = {{ preBasis_->tetrahedronOffset_ + preBasis_->dofsPerTetrahedron*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isHexahedron())
                  {
                    *it = {{ preBasis_->hexahedronOffset_ + preBasis_->dofsPerHexahedron*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isPrism())
                  {
                    *it = {{ preBasis_->prismOffset_ + preBasis_->dofsPerPrism*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isPyramid())
                  {
                    *it = {{ preBasis_->pyramidOffset_ + preBasis_->dofsPerPyramid*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
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

  auto order() const
  {
    return (k>=0) ? k : preBasis_->order();
  }

  const PreBasis* preBasis_;

  const Node* node_;
};



namespace BasisFactory {

namespace Imp {

template<int k>
class LagrangePreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  // \brief Constructor for factory with compile-time order
  LagrangePreBasisFactory()
  {}

  // \brief Constructor for factory with run-time order (template argument k is disregarded)
  LagrangePreBasisFactory(unsigned int order)
  : order_(order)
  {}

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return (k>=0)
      ? LagrangePreBasis<GridView, k, MultiIndex>(gridView)
      : LagrangePreBasis<GridView, k, MultiIndex>(gridView, order_);
  }

private:
  unsigned int order_;
};

} // end namespace BasisFactory::Imp



/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions; -1 means 'order determined at run-time'
 */
template<int k=-1>
auto lagrange()
{
  return Imp::LagrangePreBasisFactory<k>();
}

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis with a run-time order
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
auto lagrange(int order)
{
  return Imp::LagrangePreBasisFactory<-1>(order);
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
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis; -1 means 'order determined at run-time'
 */
template<typename GV, int k=-1>
using LagrangeBasis = DefaultGlobalBasis<LagrangePreBasis<GV, k, FlatMultiIndex<std::size_t>> >;





} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
