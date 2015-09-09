// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALINDEXSET_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALINDEXSET_HH



//template<typename GV, int k, class MI, class ST>
template<class LV, class NIS>
class DefaultLocalIndexSet
{
public:
  using LocalView = LV;
  using NodeIndexSet = NIS;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = typename NodeIndexSet::MultiIndex;
  using size_type = typename NodeIndexSet::size_type;


//  template<class NI>
//  DefaultLocalIndexSet(NI&& nodeIndexSet) :
//    nodeIndexSet_(std::forward<NI>(nodeIndexSet))
//  {}

  DefaultLocalIndexSet(const NodeIndexSet& nodeIndexSet) :
    nodeIndexSet_(nodeIndexSet)
  {}

  DefaultLocalIndexSet(NodeIndexSet&& nodeIndexSet) :
    nodeIndexSet_(nodeIndexSet)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const LocalView& localView)
  {
    localView_ = &localView;
    nodeIndexSet_.bind(localView_->tree());
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    localView_ = nullptr;
    nodeIndexSet_.unbind();
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return nodeIndexSet_.size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    return nodeIndexSet_.index(i);
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

protected:

  const LocalView* localView_;

  NodeIndexSet nodeIndexSet_;
};



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALINDEXSET_HH
