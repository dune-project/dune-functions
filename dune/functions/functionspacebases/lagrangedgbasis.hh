// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>




namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   LagrangeDGPreBasis
//   LagrangeDGNodeIndexSet
//   LagrangeDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k>
using LagrangeDGNode = LagrangeNode<GV, k>;

template<typename GV, int k, class MI>
class LagrangeDGNodeIndexSet;


template<typename GV, int k, class MI>
class LagrangeDGPreBasis
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


  using Node = LagrangeDGNode<GV, k>;

  using IndexSet = LagrangeDGNodeIndexSet<GV, k, MI>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  LagrangeDGPreBasis(const GridView& gv) :
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
        quadrilateralOffset_ = dofsPerTriangle * gridView_.size(Dune::GeometryTypes::triangle);
        break;
      }
      case 3:
      {
        prismOffset_         = dofsPerTetrahedron * gridView_.size(Dune::GeometryTypes::tetrahedron);

        hexahedronOffset_    = prismOffset_         +   dofsPerPrism * gridView_.size(Dune::GeometryTypes::prism);

        pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * gridView_.size(Dune::GeometryTypes::hexahedron);
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

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{};
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
    switch (dim)
    {
      case 1:
        return dofsPerEdge*gridView_.size(0);
      case 2:
      {
        return dofsPerTriangle*gridView_.size(Dune::GeometryTypes::triangle) + dofsPerQuad*gridView_.size(Dune::GeometryTypes::quadrilateral);
      }
      case 3:
      {
        return dofsPerTetrahedron*gridView_.size(Dune::GeometryTypes::tetrahedron)
             + dofsPerPyramid*gridView_.size(Dune::GeometryTypes::pyramid)
             + dofsPerPrism*gridView_.size(Dune::GeometryTypes::prism)
             + dofsPerHexahedron*gridView_.size(Dune::GeometryTypes::hexahedron);
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
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
    return StaticPower<(k+1),GV::dimension>::power;
  }

//protected:
  GridView gridView_;

  size_t quadrilateralOffset_;
  size_t pyramidOffset_;
  size_t prismOffset_;
  size_t hexahedronOffset_;
};



template<typename GV, int k, class MI>
class LagrangeDGNodeIndexSet
{
  // Cannot be an enum -- otherwise the switch statement below produces compiler warnings
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = LagrangeDGPreBasis<GV, k, MI>;

  using Node = LagrangeDGNode<GV, k>;

  LagrangeDGNodeIndexSet(const PreBasis& preBasis) :
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

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(It it) const
  {
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    for (size_type i = 0, end = size() ; i < end ; ++i, ++it)
      {
        switch (dim)
          {
          case 1:
            {
              *it = {preBasis_->dofsPerEdge*gridIndexSet.subIndex(element,0,0) + i};
              continue;
            }
          case 2:
            {
              if (element.type().isTriangle())
                {
                  *it = {preBasis_->dofsPerTriangle*gridIndexSet.subIndex(element,0,0) + i};
                  continue;
                }
              else if (element.type().isQuadrilateral())
                {
                  *it = { preBasis_->quadrilateralOffset_ + preBasis_->dofsPerQuad*gridIndexSet.subIndex(element,0,0) + i};
                  continue;
                }
              else
                DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
            }
          case 3:
            {
              if (element.type().isTetrahedron())
                {
                  *it = {preBasis_->dofsPerTetrahedron*gridIndexSet.subIndex(element,0,0) + i};
                  continue;
                }
              else if (element.type().isPrism())
                {
                  *it = { preBasis_->prismOffset_ + preBasis_->dofsPerPrism*gridIndexSet.subIndex(element,0,0) + i};
                  continue;
                }
              else if (element.type().isHexahedron())
                {
                  *it = { preBasis_->hexahedronOffset_ + preBasis_->dofsPerHexahedron*gridIndexSet.subIndex(element,0,0) + i};
                  continue;
                }
              else if (element.type().isPyramid())
                {
                  *it = { preBasis_->pyramidOffset_ + preBasis_->dofsPerPyramid*gridIndexSet.subIndex(element,0,0) + i};
                  continue;
                }
              else
                DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedrons, prisms, hexahedrons or pyramids");
            }
          }
        DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
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
class LagrangeDGPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return LagrangeDGPreBasis<GridView, k, MultiIndex>(gridView);
  }

};

} // end namespace BasisFactory::Imp



/**
 * \brief Create a pre-basis factory that can create a LagrangeDG pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions
 */
template<std::size_t k>
auto lagrangeDG()
{
  return Imp::LagrangeDGPreBasisFactory<k>();
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
 */
template<typename GV, int k>
using LagrangeDGBasis = DefaultGlobalBasis<LagrangeDGPreBasis<GV, k, FlatMultiIndex<std::size_t>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
