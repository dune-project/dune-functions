// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALINDEXSET_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALINDEXSET_HH


#include <tuple>

#include <dune/functions/functionspacebases/defaultlocalindexset.hh>



namespace Dune {
namespace Functions {



template<class LV, class NF>
class DefaultGlobalIndexSet
{
  using RootTreePath = TypeTree::HybridTreePath<>;

public:
  using LocalView = LV;
  using NodeFactory = NF;

  using NodeIndexSet = typename NodeFactory::template IndexSet<RootTreePath>;
  using SizePrefix = typename NodeFactory::SizePrefix;
  using LocalIndexSet = DefaultLocalIndexSet<LocalView, NodeIndexSet>;

  using size_type = typename NodeFactory::size_type;

  DefaultGlobalIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory)
  {}

  size_type dimension() const
  {
    return nodeFactory_->dimension();
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return nodeFactory_->size();
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return nodeFactory_->size(prefix);
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(nodeFactory_->template indexSet<RootTreePath>());
  }

private:

  const NodeFactory* nodeFactory_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALINDEXSET_HH
