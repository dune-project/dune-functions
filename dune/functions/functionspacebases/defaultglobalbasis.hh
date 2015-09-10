// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH


#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>
#include <dune/functions/functionspacebases/defaultglobalindexset.hh>
#include <dune/functions/functionspacebases/defaultlocalview.hh>



namespace Dune {
namespace Functions {



template<class NF>
class DefaultGlobalBasis
: public GridViewFunctionSpaceBasis<typename NF::GridView,
                                    DefaultLocalView<DefaultGlobalBasis<NF>>,
                                    DefaultGlobalIndexSet<DefaultLocalView<DefaultGlobalBasis<NF>>, NF>,
                                    typename NF::MultiIndex >
{
  using RootTreePath = std::tuple<>;

public:

  using NodeFactory = NF;

  //! The grid view that the FE space is defined on
  using GridView = typename NodeFactory::GridView;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = typename NodeFactory::MultiIndex;

  using size_type = typename NodeFactory::size_type;

  //! Type of the local view on the restriction of the basis to a single element
  using LocalView = DefaultLocalView<DefaultGlobalBasis<NodeFactory>>;

  using NodeIndexSet = typename NodeFactory::template IndexSet<RootTreePath>;
  using SizePrefix = typename NodeFactory::SizePrefix;
  using LocalIndexSet = DefaultLocalIndexSet<LocalView, NodeIndexSet>;
  using GlobalIndexSet = DefaultGlobalIndexSet<LocalView, NodeFactory>;


  /** \brief Constructor for a given grid view object */
  template<class... T>
  DefaultGlobalBasis(T&&... t) :
    nodeFactory_(std::forward<T>(t)...)
  {
    nodeFactory_.initializeIndices();
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return nodeFactory_.gridView();
  }

  const NodeFactory& nodeFactory() const
  {
    return nodeFactory_;
  }

  /**
   * \todo This method has been added to the interface without prior discussion.
   */
  size_type dimension() const
  {
    return nodeFactory_.dimension();
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return nodeFactory_.size();
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return nodeFactory_.size(prefix);
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(*this);
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(nodeFactory_.template indexSet<RootTreePath>());
  }

  GlobalIndexSet indexSet() const
  {
    return GlobalIndexSet(nodeFactory_);
  }

protected:
  NodeFactory nodeFactory_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH
