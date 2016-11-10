// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkNodeFactory
//   PQkNodeIndexSet
//   PQkNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename TP>
class PQkNode;

template<typename GV, int k, class MI, class TP>
class PQkNodeIndexSet;

template<typename GV, int k, class MI>
class PQkNodeFactory;



/**
 * \brief A factory for PQ-lagrange bases with given order
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
class PQkNodeFactory
{
  static const int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

private:

  template<typename, int, class, class>
  friend class PQkNodeIndexSet;

  // Precompute the number of dofs per entity type
  const static size_type dofsPerVertex =
      k == 0 ? (dim == 0 ? 1 : 0) : 1;
  const static size_type dofsPerEdge =
      k == 0 ? (dim == 1 ? 1 : 0) : k-1;
  const static size_type dofsPerTriangle =
      k == 0 ? (dim == 2 ? 1 : 0) : (k-1)*(k-2)/2;
  const static size_type dofsPerQuad =
      k == 0 ? (dim == 2 ? 1 : 0) : (k-1)*(k-1);
  const static size_type dofsPerTetrahedron =
      k == 0 ? (dim == 3 ? 1 : 0) : (k-3)*(k-2)*(k-1)/6;
  const static size_type dofsPerPrism =
      k == 0 ? (dim == 3 ? 1 : 0) : (k-1)*(k-1)*(k-2)/2;
  const static size_type dofsPerHexahedron =
      k == 0 ? (dim == 3 ? 1 : 0) : (k-1)*(k-1)*(k-1);
  const static size_type dofsPerPyramid =
      k == 0 ? (dim == 3 ? 1 : 0) : (k-2)*(k-1)*(2*k-3)/6;

public:

  //! Template mapping root tree path to type of created tree node
  template<class TP>
  using Node = PQkNode<GV, k, TP>;

  //! Template mapping root tree path to type of created tree node index set
  template<class TP>
  using IndexSet = PQkNodeIndexSet<GV, k, MI, TP>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object
  PQkNodeFactory(const GridView& gv) :
    gridView_(gv)
  {}

  //! Initialize the global indices
  void initializeIndices()
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_          + dofsPerVertex * ((size_type)gridView_.size(dim));
    triangleOffset_      = edgeOffset_            + dofsPerEdge * ((size_type) gridView_.size(dim-1));

    GeometryType triangle;
    triangle.makeTriangle();
    quadrilateralOffset_ = triangleOffset_        + dofsPerTriangle * ((size_type)gridView_.size(triangle));

    Dune::GeometryType quadrilateral;
    quadrilateral.makeQuadrilateral();
    if (dim==3) {
      tetrahedronOffset_   = quadrilateralOffset_ + dofsPerQuad * ((size_type)gridView_.size(quadrilateral));

      GeometryType tetrahedron;
      tetrahedron.makeSimplex(3);
      prismOffset_         = tetrahedronOffset_   +   dofsPerTetrahedron * ((size_type)gridView_.size(tetrahedron));

      GeometryType prism;
      prism.makePrism();
      hexahedronOffset_    = prismOffset_         +   dofsPerPrism * ((size_type)gridView_.size(prism));

      GeometryType hexahedron;
      hexahedron.makeCube(3);
      pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * ((size_type)gridView_.size(hexahedron));
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
   * \brief Create tree node with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * By passing a non-trivial root tree path this can be used
   * to create a node suitable for being placed in a tree at
   * the position specified by the root tree path.
   */
  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp};
  }

  /**
   * \brief Create tree node index set with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * Create an index set suitable for the tree node obtained
   * by node(tp).
   */
  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
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
        GeometryType triangle, quad;
        triangle.makeTriangle();
        quad.makeQuadrilateral();

        return dofsPerVertex * ((size_type)gridView_.size(dim))
          + dofsPerEdge * ((size_type)gridView_.size(dim-1))
          + dofsPerTriangle * ((size_type)gridView_.size(triangle))
          + dofsPerQuad * ((size_type)gridView_.size(quad));
      }
      case 3:
      {
        GeometryType triangle, quad, tetrahedron, pyramid, prism, hexahedron;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        tetrahedron.makeTetrahedron();
        pyramid.makePyramid();
        prism.makePrism();
        hexahedron.makeCube(3);
        return dofsPerVertex * ((size_type)gridView_.size(dim))
          + dofsPerEdge * ((size_type)gridView_.size(dim-1))
          + dofsPerTriangle * ((size_type)gridView_.size(triangle))
          + dofsPerQuad * ((size_type)gridView_.size(quad))
          + dofsPerTetrahedron * ((size_type)gridView_.size(tetrahedron))
          + dofsPerPyramid * ((size_type)gridView_.size(pyramid))
          + dofsPerPrism * ((size_type)gridView_.size(prism))
          + dofsPerHexahedron * ((size_type)gridView_.size(hexahedron));
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
    return StaticPower<(k+1),GV::dimension>::power;
  }

protected:
  GridView gridView_;

  size_type vertexOffset_;
  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type quadrilateralOffset_;
  size_type tetrahedronOffset_;
  size_type pyramidOffset_;
  size_type prismOffset_;
  size_type hexahedronOffset_;

};



template<typename GV, int k, typename TP>
class PQkNode :
  public LeafBasisNode<std::size_t, TP>
{
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<(k+1),GV::dimension>::power;

  using Base = LeafBasisNode<std::size_t,TP>;
  using FiniteElementCache = typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, k>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  PQkNode(const TreePath& treePath) :
    Base(treePath),
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



template<typename GV, int k, class MI, class TP>
class PQkNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PQkNodeFactory<GV, k, MI>;

  using Node = typename NodeFactory::template Node<TP>;

  PQkNodeIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory),
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
  MultiIndex index(size_type i) const
  {
    assert(node_ != nullptr);
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    // The dimension of the entity that the current dof is related to
    auto dofDim = dim - localKey.codim();

    if (dofDim==0) {  // vertex dof
      return {{ (size_type)(gridIndexSet.subIndex(element,localKey.subEntity(),dim)) }};
    }

    if (dofDim==1)
    {  // edge dof
      if (dim==1)  // element dof -- any local numbering is fine
        return {{ nodeFactory_->edgeOffset_
            + nodeFactory_->dofsPerEdge * ((size_type)gridIndexSet.subIndex(element,0,0))
            + localKey.index() }};
      else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        // we have to reverse the numbering if the local triangle edge is
        // not aligned with the global edge
        auto v0 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
        auto v1 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
        bool flip = (v0 > v1);
        return {{ (flip)
              ? nodeFactory_->edgeOffset_
                + nodeFactory_->dofsPerEdge*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                + (nodeFactory_->dofsPerEdge-1)-localKey.index()
              : nodeFactory_->edgeOffset_
                + nodeFactory_->dofsPerEdge*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                + localKey.index() }};
      }
    }

    if (dofDim==2)
    {
      if (dim==2)   // element dof -- any local numbering is fine
      {
        if (element.type().isTriangle())
        {
          const int interiorLagrangeNodesPerTriangle = (k-1)*(k-2)/2;
          return {{ nodeFactory_->triangleOffset_ + interiorLagrangeNodesPerTriangle*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
        }
        else if (element.type().isQuadrilateral())
        {
          const int interiorLagrangeNodesPerQuadrilateral = (k-1)*(k-1);
          return {{ nodeFactory_->quadrilateralOffset_ + interiorLagrangeNodesPerQuadrilateral*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
        }
        else
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
      } else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        if (k>3)
          DUNE_THROW(Dune::NotImplemented, "PQkNodalBasis for 3D grids is only implemented if k<=3");

        if (k==3 and !refElement.type(localKey.subEntity(), localKey.codim()).isTriangle())
          DUNE_THROW(Dune::NotImplemented, "PQkNodalBasis for 3D grids with k==3 is only implemented if the grid is a simplex grid");

        return {{ nodeFactory_->triangleOffset_ + ((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())) }};
      }
    }

    if (dofDim==3)
    {
      if (dim==3)   // element dof -- any local numbering is fine
      {
        if (element.type().isTetrahedron())
          return {{ nodeFactory_->tetrahedronOffset_ + NodeFactory::dofsPerTetrahedron*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
        else if (element.type().isHexahedron())
          return {{ nodeFactory_->hexahedronOffset_ + NodeFactory::dofsPerHexahedron*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
        else if (element.type().isPrism())
          return {{ nodeFactory_->prismOffset_ + NodeFactory::dofsPerPrism*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
        else if (element.type().isPyramid())
          return {{ nodeFactory_->pyramidOffset_ + NodeFactory::dofsPerPyramid*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
        else
          DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedra, hexahedra, prisms, or pyramids");
      } else
        DUNE_THROW(Dune::NotImplemented, "Grids of dimension larger than 3 are no supported");
    }
    DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the PQkNodalBasis");
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t k>
struct PQkNodeFactoryBuilder
{
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView>
  auto build(const GridView& gridView)
    -> PQkNodeFactory<GridView, k, MultiIndex>
  {
    return {gridView};
  }
};

} // end namespace BasisBuilder::Imp



/**
 * \brief Create a factory builder that can build a PQkNodeFactory
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of ansatz functions
 */
template<std::size_t k>
Imp::PQkNodeFactoryBuilder<k> pq()
{
  return{};
}

} // end namespace BasisBuilder



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
 * of PQkNodeFactory.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using PQkNodalBasis = DefaultGlobalBasis<PQkNodeFactory<GV, k, FlatMultiIndex<std::size_t>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
