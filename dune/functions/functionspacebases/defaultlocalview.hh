// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH


#include <tuple>



namespace Dune {
namespace Functions {



/** \brief The restriction of a finite element basis to a single element */
template<class GB>
class DefaultLocalView
{
  using RootTreePath = TypeTree::HybridTreePath<>;

public:

  //! The global FE basis that this is a view on
  using GlobalBasis = GB;

  //! The grid view the global FE basis lives on
  using GridView = typename GlobalBasis::GridView;

  //! Type of the grid element we are bound to
  using Element = typename GridView::template Codim<0>::Entity;

  //! The type used for sizes
  using size_type = typename GlobalBasis::size_type;

  //! Type used to number the global degrees of freedom
  using MultiIndex = typename GlobalBasis::MultiIndex;

  //! Tree of local finite elements / local shape function sets
  using Tree = typename GlobalBasis::NodeFactory::template Node<RootTreePath>;

  /** \brief Construct local view for a given global finite element basis */
  DefaultLocalView(const GlobalBasis& globalBasis) :
    globalBasis_(&globalBasis),
    tree_(globalBasis_->nodeFactory().node(RootTreePath()))
  {
    initializeTree(tree_);
  }

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = e;
    bindTree(tree_, element_);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    // \TODO Once we have a way to check if entities are valid, we should re-add this.
//    if (element_)
      return element_;
//    else
//      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
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

  /** \brief Total number of degrees of freedom on this element
   */
  size_type size() const
  {
    return tree_.size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   */
  size_type maxSize() const
  {
  return globalBasis_->nodeFactory().maxNodeSize();
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  const GlobalBasis* globalBasis_;
  Element element_;
  Tree tree_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH
