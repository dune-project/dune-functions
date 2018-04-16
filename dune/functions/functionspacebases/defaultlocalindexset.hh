// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALINDEXSET_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALINDEXSET_HH

#include <dune/common/std/type_traits.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/functions/functionspacebases/subspacelocalview.hh>


namespace Dune {
namespace Functions {



template<class LV, class NIS>
class DefaultLocalIndexSet
{

public:
  using LocalView = LV;
  using NodeIndexSet = NIS;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = typename LocalView::MultiIndex;
  using size_type = typename LocalView::size_type;

  DefaultLocalIndexSet() :
    localView_(nullptr)
  {}

  DefaultLocalIndexSet(const NodeIndexSet& nodeIndexSet) :
    localView_(nullptr)
  {}

  /** \brief Bind the index set to a LocalView
   */
  void bind(const LocalView& localView)
  {
    localView_ = &localView;
  }

  /** \brief Bind the index set to a SubspaceLocalView
   */
  template<class TreePath>
  void bind(const SubspaceLocalView<LocalView, TreePath>& subspaceLocalView)
  {
    bind(subspaceLocalView.rootLocalView());
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
    return localView_->size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    return localView_->index(i);
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

protected:

  const LocalView* localView_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALINDEXSET_HH
