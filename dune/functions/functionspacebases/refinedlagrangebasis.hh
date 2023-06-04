// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDLAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDLAGRANGEBASIS_HH

#include <type_traits>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/refined.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

template<typename GV, int k, typename R=double>
class RefinedLagrangeNode;

template<typename GV, int k, typename R=double>
class RefinedLagrangePreBasis;


/**
 * \brief A pre-basis for a refined Lagrange bases
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions
 * \tparam R   Range type used for shape function values
 *
 * \note This only works for simplex grids.
 */
template <typename GV, int k, typename R>
class RefinedLagrangePreBasis
{
  static const int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = RefinedLagrangeNode<GV, k, R>;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  //! Constructor for a given grid view object
  RefinedLagrangePreBasis (const GridView& gv)
    : gridView_(gv)
  {}

  //! Initialize the global indices
  void initializeIndices ()
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_   + dofsPerSimplex(0) * ((size_type)gridView_.size(dim));

    if (dim>=2)
      triangleOffset_    = edgeOffset_     + dofsPerSimplex(1) * ((size_type)gridView_.size(dim-1));

    if (dim==3)
      tetrahedronOffset_ = triangleOffset_ + dofsPerSimplex(2) * ((size_type)gridView_.size(dim-2));
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
    return Node{};
  }

  //! Same as size(prefix) with empty prefix
  size_type size () const
  {
    switch (dim)
    {
      case 1:
        return dofsPerSimplex(0) * ((size_type)gridView_.size(dim))
             + dofsPerSimplex(1) * ((size_type)gridView_.size(dim-1));
      case 2:
        return dofsPerSimplex(0) * ((size_type)gridView_.size(dim))
             + dofsPerSimplex(1) * ((size_type)gridView_.size(dim-1))
             + dofsPerSimplex(2) * ((size_type)gridView_.size(dim-2));
      case 3:
        return dofsPerSimplex(0) * ((size_type)gridView_.size(dim))
             + dofsPerSimplex(1) * ((size_type)gridView_.size(dim-1))
             + dofsPerSimplex(2) * ((size_type)gridView_.size(dim-2))
             + dofsPerSimplex(3) * ((size_type)gridView_.size(dim-3));
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Return number of possible values for next position in multi index
  template <class SizePrefix>
  size_type size (const SizePrefix& prefix) const
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
    return size();
  }

  template <typename It>
  It indices (const Node& node, It it) const
  {
    for (size_type i = 0, end = node.finiteElement().size() ; i < end ; ++it, ++i)
    {
      Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
      const auto& gridIndexSet = gridView().indexSet();
      const auto& element = node.element();

      // The dimension of the entity that the current dof is related to
      auto dofDim = dim - localKey.codim();


      if (dofDim==0)
      { //  vertex dof
        *it = {{ (size_type)(gridIndexSet.subIndex(element,localKey.subEntity(),dim)) }};
        continue;
      }

      if (dofDim==1)
      {  // edge dof
        if (dim==1)  // element dof -- any local numbering is fine
        {
          *it = {{ edgeOffset_ + dofsPerSimplex(1) * ((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
          continue;
        }
        else
        {
          const auto refElement = Dune::referenceElement<double,dim>(element.type());

          // we have to reverse the numbering if the local triangle edge is
          // not aligned with the global edge
          auto v0 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
          auto v1 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
          bool flip = (v0 > v1);
          *it = {{ (flip)
                    ? edgeOffset_
                    + dofsPerSimplex(1)*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                    + (dofsPerSimplex(1)-1)-localKey.index()
                    : edgeOffset_
                    + dofsPerSimplex(1)*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                    + localKey.index() }};
          continue;
        }
      }

      if (dofDim==2)
      {
        if (dim==2)   // element dof -- any local numbering is fine
        {
          *it = {{ triangleOffset_ + dofsPerSimplex(2)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
          continue;
        }
        else
        {
          *it = {{ triangleOffset_ + ((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())) }};
          continue;
        }
      }

      if (dofDim==3)
      {
        if (dim==3)   // element dof -- any local numbering is fine
        {
          *it = {{ tetrahedronOffset_ + dofsPerSimplex(3)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
          continue;
        }
      }
    }
    return it;
  }

  //! Polynomial order used in the local Lagrange finite-elements
  constexpr unsigned int order() const
  {
    return k;
  }

protected:
  GridView gridView_;

  //! Number of degrees of freedom assigned to a simplex (without the ones assigned to its faces!)
  size_type dofsPerSimplex (std::size_t simplexDim) const
  {
    return k == 0 ? (dim == simplexDim ? (1<<dim) : 0) : Dune::binomial(std::size_t(1),simplexDim);
  }

  size_type vertexOffset_;
  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type tetrahedronOffset_;
};



template <typename GV, int k, typename R>
class RefinedLagrangeNode
  : public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

  // refined basis only implemented for P0 and P1
  static_assert(k == 0 || k == 1);

public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = std::conditional_t<(k==0),
    Dune::RefinedP0LocalFiniteElement<typename GV::ctype,R,dim>,
    Dune::RefinedP1LocalFiniteElement<typename GV::ctype,R,dim>>;

  RefinedLagrangeNode ()
    : finiteElement_{}
    , element_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element () const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement () const
  {
    return finiteElement_;
  }

  //! Bind to element.
  void bind (const Element& e)
  {
    element_ = &e;
    this->setSize(finiteElement_.size());
  }

protected:
  const FiniteElement finiteElement_;
  const Element* element_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a  RefinedLagrange pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 */
template <int k, typename R=double>
auto refinedLagrange ()
{
  return [](const auto& gridView) {
    return RefinedLagrangePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
}

} // end namespace BasisFactory


/** \brief Nodal basis of a continuous Lagrange finite-element space on a uniformly refined simplex element
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for simplex grids.
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of RefinedLagrangeBasis.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 * \tparam R The range type of the local basis
 */
template <typename GV, int k, typename R=double>
using RefinedLagrangeBasis = DefaultGlobalBasis<RefinedLagrangePreBasis<GV,k,R> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDLAGRANGEBASIS_HH
