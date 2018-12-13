// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RANNACHERTUREKBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RANNACHERTUREKBASIS_HH

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/rannacherturek.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   RannacherTurekPreBasis
//   RannacherTurekNodeIndexSet
//   RannacherTurekNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV>
class RannacherTurekNode;

template<typename GV, class MI>
class RannacherTurekNodeIndexSet;

template<typename GV, class MI>
class RannacherTurekPreBasis;

/**
 * \brief Pre-basis for a Rannacher-Turek basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * These are Crouzeix-Raviart-elements for quadrilateral elements.
 * See
 *  Rolf Rannacher and Stefan Turek. Simple nonconforming quadrilateral Stokes
 *  element. Numerical Methods for Partial Differential Equations, 8:97–111, 1992.
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam MI  Type to be used for multi-indices
 */
template<typename GV, class MI>
class RannacherTurekPreBasis
{
  static const int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = RannacherTurekNode<GV>;

  //! Template mapping root tree path to type of created tree node index set
  using IndexSet = RannacherTurekNodeIndexSet<GV, MI>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object
  RannacherTurekPreBasis(const GridView& gv) :
    gridView_(gv)
  {}

  //! Initialize the global indices
  void initializeIndices()
  {}

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

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return (size_type)(gridView_.size(1));
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
    return 2*GV::dimension;
  }

protected:
  GridView gridView_;
};



template<typename GV>
class RannacherTurekNode :
  public LeafBasisNode
{
  static const int dim = GV::dimension;
  static const int maxSize = 2*dim;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = RannacherTurekLocalFiniteElement<typename GV::ctype, double, dim>;

  RannacherTurekNode() :
    finiteElement_(),
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
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
#ifndef NDEBUG
    if (e.type() != finiteElement_.type())
      DUNE_THROW(Dune::Exception,
                 "Rannacher-Turek elements do not exist for elements of type " << e.type());
#endif
    this->setSize(finiteElement_.size());
  }

protected:

  const FiniteElement finiteElement_;
  const Element* element_;
};



template<typename GV, class MI>
class RannacherTurekNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = RannacherTurekPreBasis<GV, MI>;

  using Node = RannacherTurekNode<GV>;

  RannacherTurekNodeIndexSet(const PreBasis& preBasis) :
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
    return (size_type)(node_->finiteElement().size());
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(It it) const
  {
    for (size_type i = 0, end = size() ; i < end ; ++i, ++it)
      {
        Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
        const auto& gridIndexSet = preBasis_->gridView().indexSet();
        const auto& element = node_->element();

        *it = {{ (size_type)(gridIndexSet.subIndex(element,localKey.subEntity(),1)) }};
      }
    return it;
  }

protected:
  const PreBasis* preBasis_;

  const Node* node_;
};



namespace BasisFactory {

namespace Imp {

template<class Dummy=void>
class RannacherTurekPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return RannacherTurekPreBasis<GridView, MultiIndex>(gridView);
  }

};

} // end namespace BasisFactory::Imp



/**
 * \brief Create a pre-basis factory that can create a Rannacher-Turek pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
template<class Dummy=void>
auto rannacherTurek()
{
  return Imp::RannacherTurekPreBasisFactory<void>();
}

} // end namespace BasisFactory




/** \brief Rannacher-Turek basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * These are Crouzeix-Raviart-elements for quadrilateral elements.
 * See
 *  Rolf Rannacher and Stefan Turek. Simple nonconforming quadrilateral Stokes
 *  element. Numerical Methods for Partial Differential Equations, 8:97–111, 1992.
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using RannacherTurekBasis = DefaultGlobalBasis<RannacherTurekPreBasis<GV, FlatMultiIndex<std::size_t> > >;

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RANNACHERTUREKBASIS_HH
