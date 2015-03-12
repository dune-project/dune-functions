// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
#include <dune/common/std/final.hh>
#else
 #ifndef DUNE_FINAL
  #define DUNE_FINAL
 #endif
#endif

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>


namespace Dune {
namespace Functions {

template<typename GV, int k>
class PQKNodalBasisLocalView;

template<typename GV, int k>
class PQKNodalBasisLeafNode;

template<typename GV, int k>
class PQKIndexSet;

template<typename GV, int k>
class PQKLocalIndexSet
{
  enum {dim = GV::dimension};

public:
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef PQKNodalBasisLocalView<GV,k> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  PQKLocalIndexSet(const PQKIndexSet<GV,k> & indexSet)
  : basisIndexSet_(indexSet)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const PQKNodalBasisLocalView<GV,k>& localView)
  {
    localView_ = &localView;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    localView_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
    return localView_->tree().finiteElement_->size();
#else
    return localView_->tree().finiteElement_->localBasis().size();
#endif
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  const MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = localView_->tree().finiteElement_->localCoefficients().localKey(i);
    const auto& gridIndexSet = basisIndexSet_.gridView_.indexSet();
    const auto& element = localView_->element();

    // The dimension of the entity that the current dof is related to
    size_t dofDim = dim - localKey.codim();

    if (dofDim==0) {  // vertex dof
      return { gridIndexSet.subIndex(element,localKey.subEntity(),dim) };
    }

    if (dofDim==1)
    {  // edge dof
      if (dim==1)   // element dof -- any local numbering is fine
        return { basisIndexSet_.edgeOffset_ + (k-1)*gridIndexSet.subIndex(element,0,0) + localKey.index() };
      else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        // we have to reverse the numbering if the local triangle edge is
        // not aligned with the global edge
        size_t v0 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
        size_t v1 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
        bool flip = (v0 > v1);
        return { (flip)
          ? basisIndexSet_.edgeOffset_ + (k-1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + (k-2)-localKey.index()
          : basisIndexSet_.edgeOffset_ + (k-1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + localKey.index() } ;
      }
    }

    if (dofDim==2)
    {
      if (dim==2)   // element dof -- any local numbering is fine
      {
        if (element.type().isTriangle())
        {
          const int interiorLagrangeNodesPerTriangle = (k-1)*(k-2)/2;
          return { basisIndexSet_.triangleOffset_ + interiorLagrangeNodesPerTriangle*gridIndexSet.subIndex(element,0,0) + localKey.index() };
        }
        else if (element.type().isQuadrilateral())
        {
          const int interiorLagrangeNodesPerQuadrilateral = (k-1)*(k-1);
          return { basisIndexSet_.quadrilateralOffset_ + interiorLagrangeNodesPerQuadrilateral*gridIndexSet.subIndex(element,0,0) + localKey.index() };
        }
        else
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
      } else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        if (! refElement.type(localKey.subEntity(), localKey.codim()).isTriangle()
            or k>3)
          DUNE_THROW(Dune::NotImplemented, "PQKNodalBasis for 3D grids is only implemented if k<=3 and if the grid is a simplex grid");

        return { basisIndexSet_.triangleOffset_ + gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) };
      }
    }
    DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the PQKNodalBasis");
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

  const PQKNodalBasisLocalView<GV,k>* localView_;

  const PQKIndexSet<GV,k> basisIndexSet_;
};

template<typename GV, int k>
class PQKIndexSet
{
  static const int dim = GV::dimension;

  // Needs the mapper
  friend class PQKLocalIndexSet<GV,k>;

  // Precompute the number of dofs per entity type
  const int dofsPerEdge        = k-1;
  const int dofsPerTriangle    = (k-1)*(k-2)/2;
  const int dofsPerQuad        = (k-1)*(k-1);
  const int dofsPerTetrahedron = (k-1)*(k-2)*(k-3)/3;
  const int dofsPerPrism       = (k-1)*(k-1)*(k-2)/2;
  const int dofsPerHexahedron  = (k-1)*(k-1)*(k-1);

public:

  typedef PQKLocalIndexSet<GV,k> LocalIndexSet;

  PQKIndexSet(const GV& gridView)
  : gridView_(gridView)
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_          + gridView_.size(dim);
    triangleOffset_      = edgeOffset_            + dofsPerEdge * gridView_.size(dim-1);

    GeometryType triangle;
    triangle.makeTriangle();
    quadrilateralOffset_ = triangleOffset_        + dofsPerTriangle * gridView_.size(triangle);

    Dune::GeometryType quadrilateral;
    quadrilateral.makeQuadrilateral();
    if (dim==3) {
      tetrahedronOffset_   = quadrilateralOffset_ + dofsPerQuad * gridView_.size(quadrilateral);

      GeometryType tetrahedron;
      tetrahedron.makeSimplex(3);
      prismOffset_         = tetrahedronOffset_   +   dofsPerTetrahedron * gridView_.size(tetrahedron);

      GeometryType prism;
      prism.makePrism();
      hexahedronOffset_    = prismOffset_         +   dofsPerPrism * gridView_.size(prism);

      GeometryType hexahedron;
      hexahedron.makeCube(3);
      pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * gridView_.size(hexahedron);
    }

  }

  std::size_t size() const
  {
    switch (dim)
    {
      case 1:
        // One for each vertex, and two for each element
        return gridView_.size(1) + dofsPerEdge*gridView_.size(0);
      case 2:
      {
        GeometryType triangle, quad;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        return gridView_.size(dim) + dofsPerEdge*gridView_.size(1)
             + gridView_.size(triangle) + dofsPerQuad*gridView_.size(quad);
      }
      case 3:
      {
        GeometryType triangle, quad, pyramid, prism, hexahedron;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        pyramid.makePyramid();
        prism.makePrism();
        hexahedron.makeCube(3);
        return gridView_.size(dim) + dofsPerEdge*gridView_.size(2)
             + gridView_.size(triangle) + dofsPerQuad*gridView_.size(quad)
             + gridView_.size(pyramid) + dofsPerPrism*gridView_.size(prism) + dofsPerHexahedron*gridView_.size(hexahedron);
      }

    }

    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(*this);
  }

private:

  size_t vertexOffset_;
  size_t edgeOffset_;
  size_t triangleOffset_;
  size_t quadrilateralOffset_;
  size_t tetrahedronOffset_;
  size_t pyramidOffset_;
  size_t prismOffset_;
  size_t hexahedronOffset_;

  const GV gridView_;
};

/** \brief Nodal basis of a scalar third-order Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - Grids must be 1d, 2d, or 3d
 * - 3d grids must be simplex grids
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
class PQKNodalBasis
: public GridViewFunctionSpaceBasis<GV,
                                    PQKNodalBasisLocalView<GV,k>,
                                    PQKIndexSet<GV,k>,
                                    std::array<std::size_t, 1> >
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef PQKNodalBasisLocalView<GV,k> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  /** \brief Constructor for a given grid view object */
  PQKNodalBasis(const GridView& gv) :
    gridView_(gv),
    indexSet_(gv)
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  PQKIndexSet<GV,k> indexSet() const
  {
    return indexSet_;
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

protected:
  const GridView gridView_;

  PQKIndexSet<GV,k> indexSet_;
};


/** \brief The restriction of a finite element basis to a single element */
template<typename GV, int k>
class PQKNodalBasisLocalView
{
public:
  /** \brief The global FE basis that this is a view on */
  typedef PQKNodalBasis<GV,k> GlobalBasis;
  typedef typename GlobalBasis::GridView GridView;

  /** \brief The type used for sizes */
  typedef typename GlobalBasis::size_type size_type;

  /** \brief Type used to number the degrees of freedom
   *
   * In the case of mixed finite elements this really can be a multi-index, but for a standard
   * P3 space this is only a single-digit multi-index, i.e., it is an integer.
   */
  typedef typename GlobalBasis::MultiIndex MultiIndex;

  /** \brief Type of the grid element we are bound to */
  typedef typename GridView::template Codim<0>::Entity Element;

  /** \brief Tree of local finite elements / local shape function sets
   *
   * In the case of a P3 space this tree consists of a single leaf only,
   * i.e., Tree is basically the type of the LocalFiniteElement
   */
  typedef PQKNodalBasisLeafNode<GV,k> Tree;

  /** \brief Construct local view for a given global finite element basis */
  PQKNodalBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(globalBasis)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = &e;
    tree_.bind(e);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   * And indeed, in the PQKNodalBasisView implementation this method does nothing.
   */
  void unbind()
  {}

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /** \brief Number of degrees of freedom on this element
   */
  size_type size() const
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
    return tree_.finiteElement_->size();
#else
    return tree_.finiteElement_->localBasis().size();
#endif
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   * The method returns k^dim, which is the number of degrees of freedom you get for cubes.
   */
  size_type maxSize() const
  {
    return StaticPower<(k+1),GV::dimension>::power;
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  const GlobalBasis* globalBasis_;
  const Element* element_;
  Tree tree_;
};


template<typename GV, int k>
class PQKNodalBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, GV::dimension, k>::FiniteElementType,
    typename PQKNodalBasis<GV,k>::size_type>
{
  typedef PQKNodalBasis<GV,k> GlobalBasis;
  static const int dim = GV::dimension;

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, k> FiniteElementCache;
  typedef typename FiniteElementCache::FiniteElementType FE;
  typedef typename GlobalBasis::size_type ST;
  typedef typename GlobalBasis::MultiIndex MI;

  typedef typename GlobalBasis::LocalView LocalView;

  friend LocalView;
  friend class PQKLocalIndexSet<GV,k>;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;

  PQKNodalBasisLeafNode(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    finiteElement_(nullptr),
    element_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const DUNE_FINAL
  {
    return *finiteElement_;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const DUNE_FINAL
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
    return finiteElement_->size();
#else
    return finiteElement_->localBasis().size();
#endif
  }

protected:

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
  }

  const GlobalBasis* globalBasis_;
  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
