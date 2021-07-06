// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDINDEXBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDINDEXBASIS_HH

#include <tuple>
#include <utility>

#include <dune/common/std/apply.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/utility.hh>

#include <dune/functions/common/staticforloop.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {
namespace Experimental {

// *****************************************************************************
// *****************************************************************************

/**
 * \brief A pre-basis transforming multi-indices
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \warning This is experimental and may be removed or
 * modified in a non-compatible way. When using this
 * functionality take care to follow the dune-functions
 * development to be aware of possible changes.
 *
 * This pre-basis wraps another pre-basis and transforms its global
 * multi-indices.
 *
 * \tparam MI  Type to be used for global multi-indices
 * \tparam RPB Raw PreBasis to be wrapped
 * \tparam T Class of the index transformation
 */
template<class MI, class RPB, class T>
class TransformedIndexPreBasis
{
  using Transformation = T;

  using This = TransformedIndexPreBasis<MI, RPB, T>;

public:

  using RawPreBasis = RPB;

  //! The grid view that the FE basis is defined on
  using GridView = typename RawPreBasis::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = typename RawPreBasis::Node;

  //! Type of created tree node index set. \deprecated
  using IndexSet = Impl::DefaultNodeIndexSet<TransformedIndexPreBasis>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, MultiIndex::max_size()+1>;

  /**
   * \brief Constructor for given child pre-basis objects
   *
   * The child pre-basis will be stored as copies
   */
  template<class RPB_R, class T_R>
  TransformedIndexPreBasis(RPB_R&& rawPreBasis, T_R&& transformation) :
    rawPreBasis_(std::forward<RPB_R>(rawPreBasis)),
    transformation_(std::forward<T_R>(transformation))
  {}

  //! Initialize the global indices
  void initializeIndices()
  {
    rawPreBasis_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return rawPreBasis_.gridView();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
    rawPreBasis_.update(gv);
  }

  /**
   * \brief Create tree node with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * By passing a non-trivial root tree path this can be used
   * to create a node suitable for being placed in a tree at
   * the position specified by the root tree path.
   */
  Node makeNode() const
  {
    return rawPreBasis_.makeNode();
  }

  /**
   * \brief Create tree node index set
   *
   * Create an index set suitable for the tree node obtained
   * by makeNode().
   * \deprecated
   */
  [[deprecated("Warning: The IndexSet typedef and the makeIndexSet method are deprecated. "\
               "As a replacement use the indices() method of the PreBasis directly.")]]
  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return size({});
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return transformation_.size(prefix, rawPreBasis_);
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return transformation_.dimension(rawPreBasis_);
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return rawPreBasis_.maxNodeSize();
  }

  const RawPreBasis& rawPreBasis() const
  {
    return rawPreBasis_;
  }

  RawPreBasis& rawPreBasis()
  {
    return rawPreBasis_;
  }

  void transformIndex(MultiIndex& multiIndex) const
  {
    transformation_.transformIndex(multiIndex, rawPreBasis_);
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    Impl::preBasisIndices(rawPreBasis(), node, it);
    for(std::size_t i=0; i<node.size(); ++i)
    {
      transformIndex(*it);
      ++it;
    }
    return it;
  }

protected:
  RawPreBasis rawPreBasis_;
  Transformation transformation_;
};



} // end namespace Experimental


namespace BasisFactory {
namespace Experimental {

namespace Imp {

template<class RawPreBasisFactory, class Transformation>
class TransformedIndexPreBasisFactory
{
public:

  static const std::size_t requiredMultiIndexSize = Transformation::maxIndexSize;

  template<class RPBF_R, class T_R>
  TransformedIndexPreBasisFactory(RPBF_R&& rawPreBasisFactory, T_R&& transformation) :
    rawPreBasisFactory_(std::forward<RPBF_R>(rawPreBasisFactory)),
    transformation_(std::forward<T_R>(transformation))
  {}

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    auto rawPreBasis = rawPreBasisFactory_.template makePreBasis<MultiIndex>(gridView);
    using RawPreBasis = std::decay_t<decltype(rawPreBasis)>;
    return Dune::Functions::Experimental::TransformedIndexPreBasis<MultiIndex, RawPreBasis, Transformation>(std::move(rawPreBasis), std::move(transformation_));
  }

private:
  RawPreBasisFactory rawPreBasisFactory_;
  Transformation transformation_;
};

} // end namespace BasisFactory::Experimental::Imp



/**
 * \brief Create a TransformedIndexPreBasisFactory
 *
 * \warning This is experimental and may be removed or
 * modified in a non-compatible way. When using this
 * functionality take care to follow the dune-functions
 * development to be aware of possible changes.
 *
 * \param preBasisFactory A PreBasisFactory creating the wrapped pre-basis
 * \param transformation The transformation object
 */
template<class RawPreBasisFactory, class Transformation>
auto transformIndices(
    RawPreBasisFactory&& preBasisFactory,
    Transformation&& transformation)
{
  return Imp::TransformedIndexPreBasisFactory<std::decay_t<RawPreBasisFactory>, std::decay_t<Transformation>>(
        std::forward<RawPreBasisFactory>(preBasisFactory),
        std::forward<Transformation>(transformation));
}



/**
 * \brief A generic implementation of a transformation
 *
 * \warning This is experimental and may be removed or
 * modified in a non-compatible way. When using this
 * functionality take care to follow the dune-functions
 * development to be aware of possible changes.
 *
 * This implements the transformation based on two callbacks: One transforms an
 * existing multi-index inplace, the other implements the size() method of the
 * pre-basis for a given prefix. Both are passed the wrapped pre-basis as second
 * argument.
 *
 * \tparam IndexTransformation Callback type for transforming multi-indices
 * \tparam SizeImplementation Callback type for implemenation of size(prefix)
 * \tparam minIS Minimal multi-index size
 * \tparam maxIS Maximal multi-index size. Notice that this has to large enough to also store the untransformed indices.
 */
template<class IndexTransformation, class SizeImplementation, std::size_t minIS, std::size_t maxIS>
class GenericIndexingTransformation
{
public:

  static constexpr std::size_t minIndexSize = minIS;
  static constexpr std::size_t maxIndexSize = maxIS;

  template<class IT_R, class SI_R>
  GenericIndexingTransformation(IT_R&& indexTransformation, SI_R&& sizeImplementation) :
    indexTransformation_(std::forward<IT_R>(indexTransformation)),
    sizeImplementation_(std::forward<SI_R>(sizeImplementation))
  {}

  template<class MultiIndex, class PreBasis>
  void transformIndex(MultiIndex& multiIndex, const PreBasis& preBasis) const
  {
    indexTransformation_(multiIndex, preBasis);
  }

  template<class Prefix, class PreBasis>
  auto size(const Prefix& prefix, const PreBasis& preBasis) const
  {
    return sizeImplementation_(prefix, preBasis);
  }

  template<class PreBasis>
  auto dimension(const PreBasis& preBasis) const
  {
    return preBasis.dimension();
  }

private:
  IndexTransformation indexTransformation_;
  SizeImplementation sizeImplementation_;
};



/**
 * \brief A generic implementation of a transformation
 *
 * \warning This is experimental and may be removed or
 * modified in a non-compatible way. When using this
 * functionality take care to follow the dune-functions
 * development to be aware of possible changes.
 *
 * This implements an index-transformation based on two callbacks: One transforms an
 * existing multi-index inplace, the other implements the size() method of the
 * pre-basis for a given prefix. Both are passed the wrapped pre-basis as second
 * argument.
 *
 * \tparam IndexTransformation Callback type for transforming multi-indices
 * \tparam SizeImplementation Callback type for implemenation of size(prefix)
 * \tparam minIS Minimal multi-index size
 * \tparam maxIS Maximal multi-index size. Notice that this has to be large enough to also store the untransformed indices.
 */
template<class IndexTransformation, class SizeImplementation, std::size_t minIndexSize, std::size_t maxIndexSize>
auto indexTransformation(IndexTransformation&& indexTransformation, SizeImplementation&& sizeImplementation, Dune::index_constant<minIndexSize>, Dune::index_constant<maxIndexSize>)
{
  return GenericIndexingTransformation<
    std::decay_t<IndexTransformation>,
    std::decay_t<SizeImplementation>,
    minIndexSize, maxIndexSize>(
        std::forward<IndexTransformation>(indexTransformation),
        std::forward<SizeImplementation>(sizeImplementation));
}


} // end namespace Experimental
} // end namespace BasisFactory
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDINDEXBASIS_HH
