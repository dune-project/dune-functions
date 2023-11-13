// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH

#include <type_traits>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismixin.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the LagrangeBasis. It contains
//
//   LagrangePreBasis
//   LagrangeNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename R=double>
class LagrangeNode;

template<typename GV, int k, typename R=double>
class LagrangePreBasis;



/**
 * \brief A pre-basis for a PQ-lagrange bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   Range type used for shape function values
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 */
template<typename GV, int k, typename R>
class LagrangePreBasis :
  public LeafPreBasisMixin< LagrangePreBasis<GV,k,R> >
{
  static const int dim = GV::dimension;
  static const bool useDynamicOrder = (k<0);

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = LagrangeNode<GV, k, R>;

  //! Constructor for a given grid view object with compile-time order
  LagrangePreBasis(const GridView& gv)
  : LagrangePreBasis(gv, std::numeric_limits<unsigned int>::max())
  {}

  //! Constructor for a given grid view object and run-time order
  LagrangePreBasis(const GridView& gv, unsigned int order) :
    gridView_(gv), order_(order)
  {
    if (!useDynamicOrder && order!=std::numeric_limits<unsigned int>::max())
      DUNE_THROW(RangeError, "Template argument k has to be -1 when supplying a run-time order!");

    for (int i=0; i<=dim; i++)
    {
      dofsPerCube_[i] = computeDofsPerCube(i);
      dofsPerSimplex_[i] = computeDofsPerSimplex(i);
    }
    dofsPerPrism_ = computeDofsPerPrism();
    dofsPerPyramid_ = computeDofsPerPyramid();
  }

  //! Initialize the global indices
  void initializeIndices()
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

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
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

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    // That cast to unsigned int is necessary because GV::dimension is an enum,
    // which is not recognized by the power method as an integer type...
    return power(order()+1, (unsigned int)GV::dimension);
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    for (size_type i = 0, end = node.finiteElement().size() ; i < end ; ++it, ++i)
    {
      Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
      const auto& gridIndexSet = gridView().indexSet();
      const auto& element = node.element();

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
              *it = {{ edgeOffset_
                       + dofsPerCube(1) * ((size_type)gridIndexSet.subIndex(element,0,0))
                       + localKey.index() }};
              continue;
            }
          else
            {
              const auto refElement
                = Dune::referenceElement<double,dim>(element.type());

              // We have to reverse the numbering if the local element edge is
              // not aligned with the global edge.
              auto v0 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
              auto v1 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
              bool flip = (v0 > v1);
              *it = {{ (flip)
                       ? edgeOffset_
                       + dofsPerCube(1)*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                       + (dofsPerCube(1)-1)-localKey.index()
                       : edgeOffset_
                       + dofsPerCube(1)*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
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
                  *it = {{ triangleOffset_ + dofsPerSimplex(2)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else if (element.type().isQuadrilateral())
                {
                  *it = {{ quadrilateralOffset_ + dofsPerCube(2)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else
                DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
            } else
            {
              const auto refElement
                = Dune::referenceElement<double,dim>(element.type());

              if (order()>3)
                DUNE_THROW(Dune::NotImplemented, "LagrangeBasis for 3D grids is only implemented if k<=3");

              if (order()==3 and !refElement.type(localKey.subEntity(), localKey.codim()).isTriangle())
                DUNE_THROW(Dune::NotImplemented, "LagrangeBasis for 3D grids with k==3 is only implemented if the grid is a simplex grid");

              *it = {{ triangleOffset_ + ((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())) }};
              continue;
            }
        }

      if (dofDim==3)
        {
          if (dim==3)   // element dof -- any local numbering is fine
            {
              if (element.type().isTetrahedron())
                {
                  *it = {{ tetrahedronOffset_ + dofsPerSimplex(3)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else if (element.type().isHexahedron())
                {
                  *it = {{ hexahedronOffset_ + dofsPerCube(3)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else if (element.type().isPrism())
                {
                  *it = {{ prismOffset_ + dofsPerPrism()*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else if (element.type().isPyramid())
                {
                  *it = {{ pyramidOffset_ + dofsPerPyramid()*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else
                DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedra, hexahedra, prisms, or pyramids");
            } else
            DUNE_THROW(Dune::NotImplemented, "Grids of dimension larger than 3 are no supported");
        }
      DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the LagrangeBasis");
    }
    return it;
  }

  //! Polynomial order used in the local Lagrange finite-elements
  unsigned int order() const
  {
    return (useDynamicOrder) ? order_ : k;
  }

protected:
  GridView gridView_;

  // Run-time order, only valid if k<0
  unsigned int order_;

  //! Number of degrees of freedom assigned to a simplex (without the ones assigned to its faces!)
  size_type dofsPerSimplex(std::size_t simplexDim) const
  {
    return useDynamicOrder ? dofsPerSimplex_[simplexDim] : computeDofsPerSimplex(simplexDim);
  }

  //! Number of degrees of freedom assigned to a cube (without the ones assigned to its faces!)
  size_type dofsPerCube(std::size_t cubeDim) const
  {
    return useDynamicOrder ? dofsPerCube_[cubeDim] : computeDofsPerCube(cubeDim);
  }

  size_type dofsPerPrism() const
  {
    return useDynamicOrder ? dofsPerPrism_ : computeDofsPerPrism();
  }

  size_type dofsPerPyramid() const
  {
    return useDynamicOrder ? dofsPerPyramid_ : computeDofsPerPyramid();
  }

  //! Number of degrees of freedom assigned to a simplex (without the ones assigned to its faces!)
  size_type computeDofsPerSimplex(std::size_t simplexDim) const
  {
    return order() == 0 ? (dim == simplexDim ? 1 : 0) : Dune::binomial(std::size_t(order()-1),simplexDim);
  }

  //! Number of degrees of freedom assigned to a cube (without the ones assigned to its faces!)
  size_type computeDofsPerCube(std::size_t cubeDim) const
  {
    return order() == 0 ? (dim == cubeDim ? 1 : 0) : Dune::power(order()-1, cubeDim);
  }

  size_type computeDofsPerPrism() const
  {
    return order() == 0 ? (dim == 3 ? 1 : 0) : (order()-1)*(order()-1)*(order()-2)/2;
  }

  size_type computeDofsPerPyramid() const
  {
    return order() == 0 ? (dim == 3 ? 1 : 0) : (order()-2)*(order()-1)*(2*order()-3)/6;
  }

  // When the order is given at run-time, the following numbers are pre-computed:
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



template<typename GV, int k, typename R>
class LagrangeNode :
  public LeafBasisNode
{
  // Stores LocalFiniteElement implementations with run-time order as a function of GeometryType
  template<typename Domain, typename Range, int dim>
  class LagrangeRunTimeLFECache
  {
  public:
    using FiniteElementType = LagrangeLocalFiniteElement<EquidistantPointSet,dim,Domain,Range>;

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

  static constexpr int dim = GV::dimension;
  static constexpr bool useDynamicOrder = (k<0);

  using FiniteElementCache = std::conditional_t<(useDynamicOrder),
                                                       LagrangeRunTimeLFECache<typename GV::ctype, R, dim>,
                                                       LagrangeLocalFiniteElementCache<typename GV::ctype, R, dim, std::max(k,0)>
                                                      >;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  //! Constructor without order (uses the compile-time value)
  LagrangeNode() :
    finiteElement_(nullptr),
    element_(nullptr)
  {}

  //! Constructor with a run-time order
  LagrangeNode(unsigned int order) :
    order_(order),
    finiteElement_(nullptr),
    element_(nullptr)
  {
    // Only the cache for the run-time-order case (i.e., k<0), has the 'order_' member
    if constexpr (useDynamicOrder)
      cache_.order_ = order;
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
    return (useDynamicOrder) ? order_ : k;
  }

  // Run-time order, only valid if k<0
  unsigned int order_;

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   The range type of the local basis
 */
template<std::size_t k, typename R=double>
auto lagrange()
{
  return [](const auto& gridView) {
    return LagrangePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
}

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis with a run-time order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 */
template<typename R=double>
auto lagrange(int order)
{
  return [=](const auto& gridView) {
    return LagrangePreBasis<std::decay_t<decltype(gridView)>, -1, R>(gridView, order);
  };
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
template<typename GV, int k=-1, typename R=double>
using LagrangeBasis = DefaultGlobalBasis<LagrangePreBasis<GV, k, R> >;





} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
