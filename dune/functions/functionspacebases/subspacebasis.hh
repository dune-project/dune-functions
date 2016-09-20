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
  using LocalIndexSet = typename RootBasis::LocalIndexSet;


  /** \brief Constructor for a given grid view object */
  SubspaceBasis(const RootBasis& rootBasis, const PrefixPath& prefixPath) :
    rootBasis_(&rootBasis),
    prefixPath_(prefixPath)
  {
//    static_assert(models<Concept::NodeFactory<GridView>, NodeFactory>(), "Type passed to DefaultGlobalBasis does not model the NodeFactory concept.");
  }

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

  LocalIndexSet localIndexSet() const
  {
    return rootBasis_->localIndexSet();
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


template<class RootBasis, class... PrefixTreeIndices>
auto subspaceBasis(const RootBasis& rootBasis, const TypeTree::HybridTreePath<PrefixTreeIndices...>& prefixPath)
{
  using PrefixPath = TypeTree::HybridTreePath<PrefixTreeIndices...>;
  return SubspaceBasis<RootBasis, PrefixPath>{rootBasis, prefixPath};
}

template<class RootBasis, class... PrefixTreeIndices>
auto subspaceBasis(const RootBasis& rootBasis, const PrefixTreeIndices&... prefixTreeIndices)
{
  return subspaceBasis(rootBasis, TypeTree::hybridTreePath(prefixTreeIndices...));
}



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH
