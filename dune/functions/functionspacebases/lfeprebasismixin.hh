// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LFEPREBASISMIXIN_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LFEPREBASISMIXIN_HH

#include <cassert>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune::Functions {

/**
 * \brief A pre-basis mixin class parametrized with a local finite-element and a DOF layout.
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * This mixin class allows for simple construction of leaf pre-bases that are based on
 * a local finite-element and a DOF layout only. Examples are the refined Lagrange pre-bases,
 * or a hierarchical Lagrange pre-basis. Note that the layout is currently not capable of
 * describing a reordering of local DOFs if there are multiple assigned to a grid entity.
 * Thus higher-order continuous finite-elements are currently not possible to describe by
 * this mixin class. Note also that this mixin fixes the local finite-element type and thus
 * cannot handle mixed GeometryTypes.
 *
 * \b Example
 * \code{.cpp}
   template <class GV, class R = double>
   class RefinedP0PreBasis :
      public LFEPreBasisMixin<GV, RefinedP0LocalFiniteElement<typename GV::ctype,R,GV::dimension>>
   {
     using LFE = RefinedP0LocalFiniteElement<typename GV::ctype,R,GV::dimension>;
     using Base = LFEPreBasisMixin<GV, LFE>;
     static const int dim = GV::dimension;
   public:
     RefinedP0PreBasis (const GV& gv) :
       Base(gv, [](GeometryType gt, int) { return (gt.dim()==dim) ? (1<<dim) : 0; })
     {}
   };
 * \endcode
 *
 * \tparam GV   The grid view that the FE basis is defined on
 * \tparam LFE  The local finite-element type
 */
template <class GV, class LFE>
class LFEPreBasisMixin :
  public LeafPreBasisMapperMixin< GV >
{
  using Base = LeafPreBasisMapperMixin< GV >;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type of the tree node
  class Node;

  /**
   * \brief Constructor for a given grid view object and layout.
   *
   * Requires that the local-finite element is default constructible.
   */
  template <class LFE_ = LFE,
    std::enable_if_t<std::is_default_constructible_v<LFE_>, int> = 0>
  LFEPreBasisMixin (const GridView& gv, MCMGLayout layout)
    : Base(gv, layout)
    , lfe_{}
  {}

  /**
   * \brief Constructor for a given grid view object, local finite-element and layout.
   *
   * Requires that the local-finite element is copyable or movable.
   */
  template <class LFE_>
  LFEPreBasisMixin (const GridView& gv, LFE_&& lfe, MCMGLayout layout)
    : Base(gv, layout)
    , lfe_(std::forward<LFE_>(lfe))
  {}

  //! Create tree node
  Node makeNode () const
  {
    return Node(lfe_);
  }

private:
  LFE lfe_;
};

// deduction guide
template <class GV, class LFE>
LFEPreBasisMixin(const GV&, const LFE&, MCMGLayout)
  -> LFEPreBasisMixin<GV,LFE>;



/**
 * \brief Leaf basis node that encapsulates a local finite-element
 * given from the LFEPreBasisMixin of type `LFE`.
 *
 * The Node implements the `LEafBasisNode` interface. Its stores a
 * pointer to the local finite-element given by the `LFEPreBasisMixin`.
 * Thus, the lifetime of the pre-basis must be greater than the lifetime
 * of this node.
 **/
template <class GV, class LFE>
class LFEPreBasisMixin<GV,LFE>::Node
  : public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = LFE;

  //! Constructor; stores a pointer to the passed local finite-element `lfe`.
  explicit Node (const LFE& lfe)
    : lfe_{&lfe}
    , element_{nullptr}
  {}

  //! Return current element; might raise an error if unbound
  const Element& element () const
  {
    assert(!!element_);
    return *element_;
  }

  /**
   * \brief Return the LocalFiniteElement for the element we are bound to; might raise an error if unbound.
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement () const
  {
    assert(!!lfe_);
    return *lfe_;
  }

  //! Bind to element. Stores a pointer to the passed element reference.
  void bind (const Element& e)
  {
    element_ = &e;
    this->setSize(lfe_->size());
  }

protected:
  const FiniteElement* lfe_;
  const Element* element_;
};


} // end namespace Dune::Functions


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LFEPREBASISMIXIN_HH
