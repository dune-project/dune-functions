// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PERIODICBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PERIODICBASIS_HH

#include <utility>
#include <type_traits>
#include <limits>
#include <set>
#include <vector>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/transformedindexbasis.hh>


namespace Dune::Functions {

namespace BasisFactory {

// The PeriodicBasis class is in the Experimental namespace because we are
// not completely sure yet whether we like it.  We reserve the right to
// modify it without advance warning.  Use at your own risk!

namespace Experimental {


/**
 * \brief Container storing identified indices for a periodic basis
 *
 * This class is intended to be passed to the BasisFactory::periodic()
 * function.
 * The class stores a set of index pairs which whould be
 * identified in order to construct a basis with periodic functions.
 */
class PeriodicIndexSet
{
  using IndexPairSet = std::set<std::pair<std::size_t,std::size_t>>;
public:

  /**
   * \brief Insert a pair of indices
   *
   * The two bases functions associated to the provided
   * indices will be combined into a single basis function
   * by associating them to a shared global index.
   */
  void unifyIndexPair(std::size_t a, std::size_t b)
  {
    if (a>b)
      std::swap(a,b);
    if (a==b)
      return;
    indexPairSet_.insert(std::make_pair(a,b));
  }

  const auto& indexPairSet() const
  {
    return indexPairSet_;
  }

private:
  IndexPairSet indexPairSet_;
};



namespace Impl {

// An index transformation for a TransformedIndexPreBasis
// impelementing periodic functions by merging indices.
// Currently only flat indices are supported.
class PeriodicIndexingTransformation
{
public:

  static constexpr std::size_t minIndexSize = 1;
  static constexpr std::size_t maxIndexSize = 1;

  template<class RawPreBasis, class IndexPairSet>
  PeriodicIndexingTransformation(const RawPreBasis& rawPreBasis, const IndexPairSet& indexPairSet)
  {
    std::size_t invalid = {std::numeric_limits<std::size_t>::max()};
    mappedIdx_.resize(rawPreBasis.size(), invalid);
    numIndices_ = 0;
    std::size_t i = 0;
    for(const auto& [a, b] : indexPairSet)
    {
      for(; i<=a; ++i)
        if (mappedIdx_[i] == invalid)
          mappedIdx_[i] = numIndices_++;
      mappedIdx_[b] = mappedIdx_[a];
    }
    for(; i<rawPreBasis.size(); ++i)
      if (mappedIdx_[i] == invalid)
        mappedIdx_[i] = numIndices_++;
  }

  template<class MultiIndex, class PreBasis>
  void transformIndex(MultiIndex& multiIndex, const PreBasis& preBasis) const
  {
    multiIndex = {{ mappedIdx_[multiIndex[0]] }};
  }

  template<class Prefix, class PreBasis>
  std::size_t size(const Prefix& prefix, const PreBasis& preBasis) const
  {
    if (prefix.size() == 1)
      return 0;
    return numIndices_;
  }

  template<class PreBasis>
  auto dimension(const PreBasis& preBasis) const
  {
    return numIndices_;
  }

private:
  std::vector<std::size_t> mappedIdx_;
  std::size_t numIndices_;
};



template<class RawPreBasisIndicator>
class PeriodicPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  PeriodicPreBasisFactory()
  {}

  template<class RPBI, class PIS>
  PeriodicPreBasisFactory(RPBI&& rawPreBasisIndicator, PIS&& periodicIndexSet) :
    rawPreBasisIndicator_(std::forward<RPBI>(rawPreBasisIndicator)),
    periodicIndexSet_(std::forward<PIS>(periodicIndexSet))
  {}

  template<class MultiIndex, class GridView,
    std::enable_if_t<models<Concept::GlobalBasis<GridView>,RawPreBasisIndicator>(), int> = 0>
  auto makePreBasis(const GridView& gridView) const
  {
    const auto& rawPreBasis = rawPreBasisIndicator_.preBasis();
    using RawPreBasis = std::decay_t<decltype(rawPreBasis)>;
    PeriodicIndexingTransformation transformation(rawPreBasis, periodicIndexSet_.indexPairSet());
    return Dune::Functions::Experimental::TransformedIndexPreBasis<MultiIndex, RawPreBasis, PeriodicIndexingTransformation>(std::move(rawPreBasis), std::move(transformation));
  }

  template<class MultiIndex, class GridView,
    std::enable_if_t<models<Concept::PreBasis<GridView>,RawPreBasisIndicator>(), int> = 0>
  auto makePreBasis(const GridView& gridView) const
  {
    const auto& rawPreBasis = rawPreBasisIndicator_;
    using RawPreBasis = std::decay_t<decltype(rawPreBasis)>;
    PeriodicIndexingTransformation transformation(rawPreBasis, periodicIndexSet_.indexPairSet());
    return Dune::Functions::Experimental::TransformedIndexPreBasis<MultiIndex, RawPreBasis, PeriodicIndexingTransformation>(std::move(rawPreBasis), std::move(transformation));
  }

  template<class MultiIndex, class GridView,
    std::enable_if_t<not models<Concept::GlobalBasis<GridView>,RawPreBasisIndicator>(), int> = 0,
    std::enable_if_t<not models<Concept::PreBasis<GridView>,RawPreBasisIndicator>(), int> = 0>
  auto makePreBasis(const GridView& gridView) const
  {
    auto rawPreBasis = rawPreBasisIndicator_.template makePreBasis<MultiIndex>(gridView);
    rawPreBasis.initializeIndices();
    using RawPreBasis = std::decay_t<decltype(rawPreBasis)>;
    PeriodicIndexingTransformation transformation(rawPreBasis, periodicIndexSet_.indexPairSet());
    return Dune::Functions::Experimental::TransformedIndexPreBasis<MultiIndex, RawPreBasis, PeriodicIndexingTransformation>(std::move(rawPreBasis), std::move(transformation));
  }

private:
  RawPreBasisIndicator rawPreBasisIndicator_;
  PeriodicIndexSet periodicIndexSet_;
};

} // end namespace BasisFactory::Impl



/**
 * \brief Create a pre-basis factory that can create a periodic pre-basis
 *
 * \param rawPreBasisIndicator Object encoding the raw non-periodic basis
 * \param periodicIndexSet A PeriodicIndexSet containing the indices to be identified
 *
 * The rawPreBasisIndicator can either be a PreBasisFactory, a PreBasis, or a GlobalBasis.
 * In the latter two cases the multi index type of those bases and the periodic
 * basis to be constructed, has to coincide. Both arguments will be copied.
 * Currently only wrapped bases with flat indices are supported.
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
template<class RawPreBasisIndicator, class PIS>
auto periodic(
    RawPreBasisIndicator&& rawPreBasisIndicator,
    PIS&& periodicIndexSet
    )
{
  return Impl::PeriodicPreBasisFactory<std::decay_t<RawPreBasisIndicator>>(
        std::forward<RawPreBasisIndicator>(rawPreBasisIndicator),
        std::forward<PIS>(periodicIndexSet));
}

} // end namespace Experimental

} // end namespace BasisFactory

} // end namespace Dune::Functions

#endif // DUNE_FUFEM_PERIODICBASIS_HH
