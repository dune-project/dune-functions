// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBSPACELOCALVIEW_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBSPACELOCALVIEW_HH


#include <tuple>

#include <dune/common/concept.hh>

#include <dune/typetree/childextraction.hh>

#include <dune/functions/functionspacebases/concepts.hh>



namespace Dune {
namespace Functions {



template<class RB, class PP>
class SubspaceBasis;



/** \brief The restriction of a finite element basis to a single element */
template<class RLV, class PP>
class SubspaceLocalView
{
  using PrefixPath = PP;

public:

  using RootLocalView = RLV;

  //! The global FE basis that this is a view on
  using GlobalBasis = SubspaceBasis<typename RootLocalView::GlobalBasis, PrefixPath>;

  //! The grid view the global FE basis lives on
  using GridView = typename GlobalBasis::GridView;

  //! Type of the grid element we are bound to
  using Element = typename GridView::template Codim<0>::Entity;

  //! The type used for sizes
  using size_type = std::size_t;

  //! Tree of local finite elements / local shape function sets
  using RootTree = typename RootLocalView::Tree;

  //! Tree of local finite elements / local shape function sets
  using Tree = typename TypeTree::ChildForTreePath<RootTree, PrefixPath>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = typename RootLocalView::MultiIndex;

  /** \brief Construct local view for a given global finite element basis */
  SubspaceLocalView(const GlobalBasis& globalBasis, const PrefixPath& /*prefixPath*/) :
    globalBasis_(&globalBasis),
    rootLocalView_(globalBasis.rootBasis().localView())
  {
//    static_assert(models<Concept::BasisTree<GridView>, Tree>(), "Tree type passed to SubspaceLocalView does not model the BasisNode concept.");
  }

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    rootLocalView_.bind(e);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    return rootLocalView_.element();
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   */
  void unbind()
  {
    rootLocalView_.unbind();
  }

  /** \brief Return if the view is bound to a grid element
   */
  bool bound() const
  {
    return rootLocalView_.bound();
  }

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return TypeTree::child(rootLocalView_.tree(), globalBasis_->prefixPath());
  }

  /** \brief Total number of degrees of freedom on this element
   */
  size_type size() const
  {
    return rootLocalView_.size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   */
  size_type maxSize() const
  {
    return rootLocalView_.maxSize();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    return rootLocalView_.index(i);
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

  const RootLocalView& rootLocalView() const
  {
    return rootLocalView_;
  }

protected:
  const GlobalBasis* globalBasis_;
  RootLocalView rootLocalView_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBSPACELOCALVIEW_HH
