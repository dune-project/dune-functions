// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBSPACEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBSPACEBASIS_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/concept.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/subspacelocalview.hh>
#include <dune/functions/functionspacebases/concepts.hh>



namespace Dune {
namespace Functions {



namespace Impl {

  template<class... Inner, class... Outer>
  auto joinTreePaths(const TypeTree::HybridTreePath<Inner...>& inner, const TypeTree::HybridTreePath<Outer...> outer)
  {
    return TypeTree::HybridTreePath<Inner..., Outer...>(std::tuple_cat(inner._data, outer._data));
  }

  template<class InnerTP, class OuterTP>
  using JoinTreePath_t = std::decay_t<decltype(joinTreePaths(std::declval<InnerTP>(), std::declval<OuterTP>()))>;

}



template<class RB, class TP>
class SubspaceBasis
{
public:

  using RootBasis = RB;

  using RootLocalView = typename RootBasis::LocalView;

  using PrefixPath = TP;

  //! The grid view that the FE space is defined on
  using GridView = typename RootBasis::GridView;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = typename RootBasis::MultiIndex;

  using size_type = std::size_t;

  //! Type of the local view on the restriction of the basis to a single element
  using LocalView = SubspaceLocalView<RootLocalView, PrefixPath>;

  using SizePrefix = typename RootBasis::SizePrefix;


  /** \brief Constructor for a given grid view object */
  SubspaceBasis(const RootBasis& rootBasis, const PrefixPath& prefixPath) :
    rootBasis_(&rootBasis),
    prefixPath_(prefixPath)
  {}

  /** \brief Constructor from another SubspaceBasis
   *
   * This will join the tree paths and use them relative to
   * rootBasis.rootBasis().
   */
  template<class RootRootBasis, class InnerTP, class OuterTP>
  SubspaceBasis(const SubspaceBasis<RootRootBasis, InnerTP>& rootBasis, const OuterTP& prefixPath) :
    SubspaceBasis(rootBasis.rootBasis(), Impl::joinTreePaths(rootBasis.prefixPath(), prefixPath))
  {}


  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return rootBasis_->gridView();
  }

  /**
   * \todo This method has been added to the interface without prior discussion.
   */
  size_type dimension() const
  {
    return rootBasis_->dimension();
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return rootBasis_->size();
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return rootBasis_->size(prefix);
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(*this, prefixPath_);
  }

  const RootBasis& rootBasis() const
  {
    return *rootBasis_;
  }

  const PrefixPath& prefixPath() const
  {
    return prefixPath_;
  }

protected:
  const RootBasis* rootBasis_;
  PrefixPath prefixPath_;
};


// CTAD guide for a non-SubspaceBasis root basis
template<class RB, class TP>
SubspaceBasis(const RB&, const TP) -> SubspaceBasis<RB, TP>;

// CTAD guide for a SubspaceBasis root basis
template<class RootRootBasis, class InnerTP, class OuterTP>
SubspaceBasis(const SubspaceBasis<RootRootBasis, InnerTP>& rootBasis, const OuterTP& prefixPath)
  -> SubspaceBasis<std::decay_t<decltype(rootBasis.rootBasis())>, Impl::JoinTreePath_t<InnerTP, OuterTP>>;



/**
 * \brief Create SubspaceBasis from a root basis and a prefixPath
 *
 * This will not return a nested SubspaceBasis if rootBasis is already
 * a SubspaceBasis. Instead it will join the tree paths and use them
 * to construct a non-nested SubspaceBasis relative to rootBasis.rootBasis().
 *
 * \param rootBasis Create a subspace basis relative to this basis
 * \param prefixPath A prefix path of the subspace within the root basis
 */
template<class RootBasis, class... PrefixTreeIndices>
auto subspaceBasis(const RootBasis& rootBasis, const TypeTree::HybridTreePath<PrefixTreeIndices...>& prefixPath)
{
  return SubspaceBasis(rootBasis, prefixPath);
}

template<class RootBasis, class... PrefixTreeIndices>
auto subspaceBasis(const RootBasis& rootBasis, const PrefixTreeIndices&... prefixTreeIndices)
{
  return subspaceBasis(rootBasis, TypeTree::hybridTreePath(prefixTreeIndices...));
}



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH
