// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/leafprebasismixin.hh>




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

template<typename GV, int k>
using LagrangeDGNode = LagrangeNode<GV, k>;

template<typename GV, int k>
class LagrangeDGPreBasis :
  public LeafPreBasisMixin< LagrangeDGPreBasis<GV,k> >
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

  size_type dimension() const
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

  size_type maxNodeSize() const
  {
    return Dune::power(k+1, int(GV::dimension));
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    const auto& gridIndexSet = gridView().indexSet();
    const auto& element = node.element();

    size_type offset = 0;
    if constexpr (dim==1)
      offset = dofsPerEdge*gridIndexSet.subIndex(element,0,0);
    else if constexpr (dim==2)
    {
      if (element.type().isTriangle())
        offset = dofsPerTriangle*gridIndexSet.subIndex(element,0,0);
      else if (element.type().isQuadrilateral())
        offset = quadrilateralOffset_ + dofsPerQuad*gridIndexSet.subIndex(element,0,0);
      else
        DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
    }
    else if constexpr (dim==3)
    {
      if (element.type().isTetrahedron())
        offset = dofsPerTetrahedron*gridIndexSet.subIndex(element,0,0);
      else if (element.type().isPrism())
        offset = prismOffset_ + dofsPerPrism*gridIndexSet.subIndex(element,0,0);
      else if (element.type().isHexahedron())
        offset = hexahedronOffset_ + dofsPerHexahedron*gridIndexSet.subIndex(element,0,0);
      else if (element.type().isPyramid())
        offset = pyramidOffset_ + dofsPerPyramid*gridIndexSet.subIndex(element,0,0);
      else
        DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedrons, prisms, hexahedrons or pyramids");
    }
    else
      DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
    for (size_type i = 0, end = node.size() ; i < end ; ++i, ++it)
      *it = {offset + i};
    return it;
  }

  //! Polynomial order used in the local Lagrange finite-elements
  unsigned int order() const
  {
    return k;
  }

protected:
  GridView gridView_;

  size_t quadrilateralOffset_;
  size_t pyramidOffset_;
  size_t prismOffset_;
  size_t hexahedronOffset_;
};



namespace BasisFactory {

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
  return [](const auto& gridView) {
    return LagrangeDGPreBasis<std::decay_t<decltype(gridView)>, k>(gridView);
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
 */
template<typename GV, int k>
using LagrangeDGBasis = DefaultGlobalBasis<LagrangeDGPreBasis<GV, k> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
