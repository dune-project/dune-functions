// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>




namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   LagrangeDGNodeFactory
//   LagrangeDGNodeIndexSet
//   LagrangeDGNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename TP>
using LagrangeDGNode = PQkNode<GV, k, TP>;

template<typename GV, int k, class MI, class TP>
class LagrangeDGNodeIndexSet;


template<typename GV, int k, class MI>
class LagrangeDGNodeFactory
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;


  // Precompute the number of dofs per entity type
  const static int dofsPerEdge        = k+1;
  const static int dofsPerTriangle    = (k+1)*(k+2)/2;
  const static int dofsPerQuad        = (k+1)*(k+1);
  const static int dofsPerTetrahedron = (k+1)*(k+2)*(k+3)/6;
  const static int dofsPerPrism       = (k+1)*(k+1)*(k+2)/2;
  const static int dofsPerHexahedron  = (k+1)*(k+1)*(k+1);
  const static int dofsPerPyramid     = (k+1)*(k+2)*(2*k+3)/6;


  template<class TP>
  using Node = LagrangeDGNode<GV, k, TP>;

  template<class TP>
  using IndexSet = LagrangeDGNodeIndexSet<GV, k, MI, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 2>;

  /** \brief Constructor for a given grid view object */
  LagrangeDGNodeFactory(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    switch (dim)
    {
      case 1:
      {
        break;
      }
      case 2:
      {
        GeometryType triangle;
        triangle.makeTriangle();
        quadrilateralOffset_ = dofsPerTriangle * gridView_.size(triangle);
        break;
      }
      case 3:
      {
        GeometryType tetrahedron;
        tetrahedron.makeSimplex(3);
        prismOffset_         = dofsPerTetrahedron * gridView_.size(tetrahedron);

        GeometryType prism;
        prism.makePrism();
        hexahedronOffset_    = prismOffset_         +   dofsPerPrism * gridView_.size(prism);

        GeometryType hexahedron;
        hexahedron.makeCube(3);
        pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * gridView_.size(hexahedron);
        break;
      }
    }
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  void update(const GridView& gv)
  {
    gridView_ = gv;
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    switch (dim)
    {
      case 1:
        return dofsPerEdge*gridView_.size(0);
      case 2:
      {
        GeometryType triangle, quad;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        return dofsPerTriangle*gridView_.size(triangle) + dofsPerQuad*gridView_.size(quad);
      }
      case 3:
      {
        GeometryType tetrahedron, pyramid, prism, hexahedron;
        tetrahedron.makeTetrahedron();
        pyramid.makePyramid();
        prism.makePrism();
        hexahedron.makeCube(3);
        return dofsPerTetrahedron*gridView_.size(tetrahedron) + dofsPerPyramid*gridView_.size(pyramid)
             + dofsPerPrism*gridView_.size(prism) + dofsPerHexahedron*gridView_.size(hexahedron);
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return size();
    assert(false);
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return StaticPower<(k+1),GV::dimension>::power;
  }

//protected:
  GridView gridView_;

  size_t quadrilateralOffset_;
  size_t pyramidOffset_;
  size_t prismOffset_;
  size_t hexahedronOffset_;
};



template<typename GV, int k, class MI, class TP>
class LagrangeDGNodeIndexSet
{
  // Cannot be an enum -- otherwise the switch statement below produces compiler warnings
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = LagrangeDGNodeFactory<GV, k, MI>;

  using Node = typename NodeFactory::template Node<TP>;

  LagrangeDGNodeIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory)
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

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    switch (dim)
    {
      case 1:
      {
        return {nodeFactory_->dofsPerEdge*gridIndexSet.subIndex(element,0,0) + i};
      }
      case 2:
      {
        if (element.type().isTriangle())
        {
          return {nodeFactory_->dofsPerTriangle*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isQuadrilateral())
        {
          return { nodeFactory_->quadrilateralOffset_ + nodeFactory_->dofsPerQuad*gridIndexSet.subIndex(element,0,0) + i};
        }
        else
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
      }
      case 3:
      {
        if (element.type().isTetrahedron())
        {
          return {nodeFactory_->dofsPerTetrahedron*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isPrism())
        {
          return { nodeFactory_->prismOffset_ + nodeFactory_->dofsPerPrism*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isHexahedron())
        {
          return { nodeFactory_->hexahedronOffset_ + nodeFactory_->dofsPerHexahedron*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isPyramid())
        {
          return { nodeFactory_->pyramidOffset_ + nodeFactory_->dofsPerPyramid*gridIndexSet.subIndex(element,0,0) + i};
        }
        else
          DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedrons, prisms, hexahedrons or pyramids");
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Lagrangean-DG finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using LagrangeDGBasis = DefaultGlobalBasis<LagrangeDGNodeFactory<GV, k, FlatMultiIndex<std::size_t>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
