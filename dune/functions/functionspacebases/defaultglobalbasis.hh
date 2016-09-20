// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/concept.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/defaultglobalindexset.hh>
#include <dune/functions/functionspacebases/defaultlocalview.hh>
#include <dune/functions/functionspacebases/concepts.hh>



namespace Dune {
namespace Functions {



template<class NF>
class DefaultGlobalBasis
{
public:

  using NodeFactory = NF;

  using PrefixPath = TypeTree::HybridTreePath<>;

  //! The grid view that the FE space is defined on
  using GridView = typename NodeFactory::GridView;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = typename NodeFactory::MultiIndex;

  using size_type = typename NodeFactory::size_type;

  //! Type of the local view on the restriction of the basis to a single element
  using LocalView = DefaultLocalView<DefaultGlobalBasis<NodeFactory>>;

  using NodeIndexSet = typename NodeFactory::template IndexSet<PrefixPath>;
  using SizePrefix = typename NodeFactory::SizePrefix;
  using LocalIndexSet = DefaultLocalIndexSet<LocalView, NodeIndexSet>;
  using GlobalIndexSet = DefaultGlobalIndexSet<LocalView, NodeFactory>;


  /** \brief Constructor for a given grid view object */
  template<class... T,
    disableCopyMove<DefaultGlobalBasis, T...> = 0,
    enableIfConstructible<NodeFactory, T...> = 0>
  DefaultGlobalBasis(T&&... t) :
    nodeFactory_(std::forward<T>(t)...),
    prefixPath_()
  {
    static_assert(models<Concept::NodeFactory<GridView>, NodeFactory>(), "Type passed to DefaultGlobalBasis does not model the NodeFactory concept.");
    nodeFactory_.initializeIndices();
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
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
    return LocalIndexSet(nodeFactory_.template indexSet<PrefixPath>());
  }

  GlobalIndexSet indexSet() const
  {
    return GlobalIndexSet(nodeFactory_);
  }

  const DefaultGlobalBasis& rootBasis() const
  {
    return *this;
  }

  const PrefixPath& prefixPath() const
  {
    return prefixPath_;
  }

protected:
  NodeFactory nodeFactory_;
  PrefixPath prefixPath_;
};



namespace BasisBuilder {

template<class GridView, class FactoryTag>
auto makeBasis(const GridView& gridView, FactoryTag&& factoryTag)
  -> DefaultGlobalBasis<decltype(factoryTag.template build<typename Dune::ReservedVector<std::size_t, FactoryTag::requiredMultiIndexSize> >(gridView))>
{
  using MultiIndex = typename Dune::ReservedVector<std::size_t, FactoryTag::requiredMultiIndexSize>;
  return {factoryTag.template build<MultiIndex>(gridView)};
}

template<class MultiIndex, class GridView, class FactoryTag>
auto makeBasis(const GridView& gridView, FactoryTag&& factoryTag)
  -> DefaultGlobalBasis<decltype(factoryTag.template build<MultiIndex>(gridView))>
{
  return {factoryTag.template build<MultiIndex>(gridView)};
}

} // end namespace BasisBuilder



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH
