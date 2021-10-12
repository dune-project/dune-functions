// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RANNACHERTUREKBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RANNACHERTUREKBASIS_HH

#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/localfunctions/common/localfiniteelementvariant.hh>
#include <dune/localfunctions/rannacherturek.hh>
#include <dune/localfunctions/crouzeixraviart.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   RannacherTurekPreBasis
//   RannacherTurekNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV>
class RannacherTurekNode;

template<typename GV>
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
 */
template<typename GV>
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

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  //! Constructor for a given grid view object
  RannacherTurekPreBasis(const GridView& gv) :
    gridView_(gv)
  {
    for(auto type : gv.indexSet().types(0))
      if (!type.isSimplex() && !type.isCube())
        DUNE_THROW(Dune::NotImplemented, "Rannacher-Turek or Crouzeix-Raviart elements are only implemented for grids with simplex or cube elements.");
  }

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

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return (size_type)(gridView_.size(1));
  }

  //! Return number of possible values for next position in multi index
  template<class SizePrefix>
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

  template<typename It>
  It indices(const Node& node, It it) const
  {
    for (size_type i = 0, end = node.size() ; i < end ; ++i, ++it)
      {
        Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
        const auto& gridIndexSet = gridView().indexSet();
        const auto& element = node.element();

        *it = {{ (size_type)(gridIndexSet.subIndex(element,localKey.subEntity(),1)) }};
      }
    return it;
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

  constexpr static bool hasFixedElementType = Capabilities::hasSingleGeometryType<typename GV::Grid>::v;

  using CubeFiniteElement    = RannacherTurekLocalFiniteElement<typename GV::ctype,double,dim>;
  using SimplexFiniteElement = CrouzeixRaviartLocalFiniteElement<typename GV::ctype,double,dim>;

  constexpr static unsigned int  topologyId = Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId;  // meaningless if hasFixedElementType is false
  constexpr static GeometryType type = GeometryType(topologyId, GV::dimension);

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = std::conditional_t<hasFixedElementType,
                                         std::conditional_t<type.isCube(),CubeFiniteElement,SimplexFiniteElement>,
                                         LocalFiniteElementVariant<CubeFiniteElement, SimplexFiniteElement> >;

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
    if constexpr (!hasFixedElementType)
      finiteElement_ = e.type().isCube() ? static_cast<FiniteElement>(CubeFiniteElement())
                                         : static_cast<FiniteElement>(SimplexFiniteElement()) ;
    this->setSize(finiteElement_.size());
  }

protected:

  FiniteElement finiteElement_;
  const Element* element_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a Rannacher-Turek pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
template<class Dummy=void>
auto rannacherTurek()
{
  return [](const auto& gridView) {
    return RannacherTurekPreBasis<std::decay_t<decltype(gridView)>>(gridView);
  };
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
using RannacherTurekBasis = DefaultGlobalBasis<RannacherTurekPreBasis<GV> >;

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RANNACHERTUREKBASIS_HH
