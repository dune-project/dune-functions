// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>

#include <dune/typetree/powernode.hh>
#include <dune/typetree/utility.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>



namespace Dune {
namespace Functions {


// *****************************************************************************
// This is the reusable part of the power bases. It contains
//
//   PowerPreBasis
//   PowerNodeIndexSet
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<class PB, class IMS>
class PowerNodeIndexSet;



/**
 * \brief A pre-basis for power bases
 *
 * This pre-basis represents a power of a given pre-basis.
 * Its node type is a PowerBasisNodes for the given subnode.
 *
 * \tparam MI  Type to be used for multi-indices
 * \tparam IMS An IndexMergingStrategy used to merge the global indices of the child factories
 * \tparam SPB  The child pre-basis
 * \tparam C   The exponent of the power node
 */
template<class MI, class IMS, class SPB, std::size_t C>
class PowerPreBasis
{
  static const std::size_t children = C;

  template<class, class>
  friend class PowerNodeIndexSet;

public:

  //! The child pre-basis
  using SubPreBasis = SPB;

  //! The grid view that the FE basis is defined on
  using GridView = typename SPB::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Strategy used to merge the global indices of the child factories
  using IndexMergingStrategy = IMS;

  using SubNode = typename SubPreBasis::Node;

  using SubIndexSet = typename SubPreBasis::IndexSet;

  //! Template mapping root tree path to type of created tree node
  using Node = PowerBasisNode<SubNode, children>;

  //! Template mapping root tree path to type of created tree node index set
  using IndexSet = PowerNodeIndexSet<PowerPreBasis, IMS>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, MultiIndex::max_size()>;

private:

  using SubMultiIndex = MI;

public:

  /**
   * \brief Constructor for given child pre-basis objects
   *
   * The child factories will be stored as copies
   */
  template<class... SFArgs,
    disableCopyMove<PowerPreBasis, SFArgs...> = 0,
    enableIfConstructible<SubPreBasis, SFArgs...> = 0>
  PowerPreBasis(SFArgs&&... sfArgs) :
    subPreBasis_(std::forward<SFArgs>(sfArgs)...)
  {
    static_assert(models<Concept::PreBasis<GridView>, SubPreBasis>(), "Subprebasis passed to PowerPreBasis does not model the PreBasis concept.");
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    subPreBasis_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return subPreBasis_.gridView();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
    subPreBasis_.update(gv);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    auto node = Node{};
    for (std::size_t i=0; i<children; ++i)
      node.setChild(i, subPreBasis_.makeNode());
    return node;
  }

  /**
   * \brief Create tree node index set
   *
   * Create an index set suitable for the tree node obtained
   * by makeNode().
   */
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
    return size(prefix, IndexMergingStrategy{});
  }

private:

  size_type size(const SizePrefix& prefix, BasisFactory::FlatInterleaved) const
  {
    // The root index size is the root index size of a single subnode
    // multiplied by the number of subnodes because we enumerate all
    // child indices in a row.
    if (prefix.size() == 0)
      return children*subPreBasis_.size({});

    // The first prefix entry refers to one of the (root index size)
    // subindex trees. Hence we have to first compute the corresponding
    // prefix entry for a single subnode subnode. The we can append
    // the other prefix entries unmodified, because the index tree
    // looks the same after the first level.
    typename SubPreBasis::SizePrefix subPrefix;
    subPrefix.push_back(prefix[0] / children);
    for(std::size_t i=1; i<prefix.size(); ++i)
      subPrefix.push_back(prefix[i]);
    return subPreBasis_.size(subPrefix);
  }

  size_type size(const SizePrefix& prefix, BasisFactory::FlatLexicographic) const
  {
    // The size at the index tree root is the size of at the index tree
    // root of a single subnode multiplied by the number of subnodes
    // because we enumerate all child indices in a row.
    if (prefix.size() == 0)
      return children*subPreBasis_.size({});

    // The first prefix entry refers to one of the (root index size)
    // subindex trees. Hence we have to first compute the corresponding
    // prefix entry for a single subnode subnode. The we can append
    // the other prefix entries unmodified, because the index tree
    // looks the same after the first level.
    typename SubPreBasis::SizePrefix subPrefix;
    subPrefix.push_back(prefix[0] % children);
    for(std::size_t i=1; i<prefix.size(); ++i)
      subPrefix.push_back(prefix[i]);
    return subPreBasis_.size(subPrefix);
  }

  size_type size(const SizePrefix& prefix, BasisFactory::BlockedLexicographic) const
  {
    if (prefix.size() == 0)
      return children;
    typename SubPreBasis::SizePrefix subPrefix;
    for(std::size_t i=1; i<prefix.size(); ++i)
      subPrefix.push_back(prefix[i]);
    return subPreBasis_.size(subPrefix);
  }

  size_type size(const SizePrefix& prefix, BasisFactory::BlockedInterleaved) const
  {
    if (prefix.size() == 0)
      return subPreBasis_.size();

    typename SubPreBasis::SizePrefix subPrefix;
    for(std::size_t i=0; i<prefix.size()-1; ++i)
      subPrefix.push_back(prefix[i]);

    size_type r = subPreBasis_.size(subPrefix);
    if (r==0)
      return 0;
    subPrefix.push_back(prefix.back());
    r = subPreBasis_.size(subPrefix);
    if (r==0)
      return children;
    return r;
  }

public:

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return subPreBasis_.dimension() * children;
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return subPreBasis_.maxNodeSize() * children;
  }

  //! Const access to the stored prebasis of the factor in the power space
  const SubPreBasis& subPreBasis() const
  {
    return subPreBasis_;
  }

  //! Mutable access to the stored prebasis of the factor in the power space
  SubPreBasis& subPreBasis()
  {
    return subPreBasis_;
  }

private:
  SubPreBasis subPreBasis_;
};



template<class PB, class IMS>
class PowerNodeIndexSet
{
public:

  using size_type = std::size_t;
  using PreBasis = PB;
  using MultiIndex = typename PreBasis::MultiIndex;
  using Node = typename PreBasis::Node;

protected:

  using IndexMergingStrategy = IMS;
  using SubIndexSet = typename PreBasis::SubPreBasis::IndexSet;
  static const std::size_t children = PreBasis::children;

public:

  PowerNodeIndexSet(const PreBasis & preBasis) :
    preBasis_(&preBasis),
    subNodeIndexSet_(preBasis_->subPreBasis().makeIndexSet())
  {}

  void bind(const Node& node)
  {
    using namespace TypeTree::Indices;
    node_ = &node;
    subNodeIndexSet_.bind(node.child(_0));
  }

  void unbind()
  {
    node_ = nullptr;
    subNodeIndexSet_.unbind();
  }

  size_type size() const
  {
    return node_->size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(It it) const
  {
    return indices(it, IndexMergingStrategy{});
  }

private:

  template<typename It>
  It indices(It multiIndices, BasisFactory::FlatInterleaved) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    // Fill indices for first child at the beginning.
    auto next = subNodeIndexSet_.indices(multiIndices);
    // Multiply first component of all indices for first child by
    // number of children to strech the index range for interleaving.
    for (std::size_t i = 0; i<subTreeSize; ++i)
      multiIndices[i][0] *= children;
    for (std::size_t child = 1; child<children; ++child)
    {
      for (std::size_t i = 0; i<subTreeSize; ++i)
      {
        // Copy indices from first child for all other children
        // and shift them by child index to interleave indices.
        //    multiIndices[child*subTreeSize+i] = multiIndices[i];
        //    multiIndices[child*subTreeSize+i][0] = multiIndices[i][0]+child;
        (*next) = multiIndices[i];
        (*next)[0] = multiIndices[i][0]+child;
        ++next;
      }
    }
    return next;
  }

  template<typename It>
  It indices(It multiIndices, BasisFactory::FlatLexicographic) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    size_type firstIndexEntrySize = preBasis_->subPreBasis().size({});
    // Fill indices for first child at the beginning.
    auto next = subNodeIndexSet_.indices(multiIndices);
    for (std::size_t child = 1; child<children; ++child)
    {
      for (std::size_t i = 0; i<subTreeSize; ++i)
      {
        // Copy indices from first child for all other children
        // and shift them by suitable offset to get lexicographic indices.
        //     multiIndices[child*subTreeSize+i] = multiIndices[i];
        //     multiIndices[child*subTreeSize+i][0] += child*firstIndexEntrySize;
        (*next) = multiIndices[i];
        (*next)[0] += child*firstIndexEntrySize;
        ++next;
      }
    }
    return next;
  }

  static void multiIndexPushFront(MultiIndex& M, size_type M0)
  {
    M.resize(M.size()+1);
    for(std::size_t i=M.size()-1; i>0; --i)
      M[i] = M[i-1];
    M[0] = M0;
  }

  template<typename It>
  It indices(It multiIndices, BasisFactory::BlockedLexicographic) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    // Fill indices for first child at the beginning.
    auto next = subNodeIndexSet_.indices(multiIndices);
    // Insert 0 before first component of all indices for first child.
    for (std::size_t i = 0; i<subTreeSize; ++i)
      multiIndexPushFront(multiIndices[i], 0);
    for (std::size_t child = 1; child<children; ++child)
    {
      for (std::size_t i = 0; i<subTreeSize; ++i)
      {
        // Copy indices from first child for all other children and overwrite
        // zero in first component as inserted above by child index.
        //     multiIndices[child*subTreeSize+i] = multiIndices[i];
        //     multiIndices[child*subTreeSize+i][0] = child;
        (*next) = multiIndices[i];
        (*next)[0] = child;
        ++next;
      }
    }
    return next;
  }

  template<typename It>
  It indices(It multiIndices, BasisFactory::BlockedInterleaved) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    // Fill indices for first child at the beginning.
    auto next = subNodeIndexSet_.indices(multiIndices);
    // Append 0 after last component of all indices for first child.
    for (std::size_t i = 0; i<subTreeSize; ++i)
      multiIndices[i].push_back(0);
    for (std::size_t child = 1; child<children; ++child)
    {
      for (std::size_t i = 0; i<subTreeSize; ++i)
      {
        // Copy indices from first child for all other children and overwrite
        // zero in last component as appended above by child index.
        (*next) = multiIndices[i];
        (*next).back() = child;
        ++next;
      }
    }
    return next;
  }

  const PreBasis* preBasis_;
  SubIndexSet subNodeIndexSet_;
  const Node* node_;
};



namespace BasisFactory {

namespace Imp {

template<std::size_t k, class IndexMergingStrategy, class ChildPreBasisFactory>
class PowerPreBasisFactory
{
  static const bool isBlocked = std::is_same<IndexMergingStrategy,BlockedLexicographic>::value or std::is_same<IndexMergingStrategy,BlockedInterleaved>::value;

  static const std::size_t maxChildIndexSize = ChildPreBasisFactory::requiredMultiIndexSize;

public:

  static const std::size_t requiredMultiIndexSize = isBlocked ? (maxChildIndexSize+1) : maxChildIndexSize;

  PowerPreBasisFactory(const ChildPreBasisFactory& childPreBasisFactory) :
    childPreBasisFactory_(childPreBasisFactory)
  {}

  PowerPreBasisFactory(ChildPreBasisFactory&& childPreBasisFactory) :
    childPreBasisFactory_(std::move(childPreBasisFactory))
  {}

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    auto childPreBasis = childPreBasisFactory_.template makePreBasis<MultiIndex>(gridView);
    using ChildPreBasis = decltype(childPreBasis);

    return PowerPreBasis<MultiIndex,  IndexMergingStrategy, ChildPreBasis, k>(std::move(childPreBasis));
  }

private:
  ChildPreBasisFactory childPreBasisFactory_;
};

} // end namespace BasisFactory::Imp



/**
 * \brief Create a pre-basis factory that can build a PowerPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam ChildPreBasisFactory Types of child pre-basis factory
 * \tparam IndexMergingStrategy An IndexMergingStrategy type
 * \param childPreBasisFactory Child pre-basis factory
 * \param ims IndexMergingStrategy to be used
 *
 * This overload can be used to explicitly supply an IndexMergingStrategy.
 */
template<std::size_t k, class ChildPreBasisFactory, class IndexMergingStrategy>
auto power(ChildPreBasisFactory&& childPreBasisFactory, const IndexMergingStrategy& ims)
{
  return Imp::PowerPreBasisFactory<k, IndexMergingStrategy, ChildPreBasisFactory>(std::forward<ChildPreBasisFactory>(childPreBasisFactory));
}

/**
 * \brief Create a factory builder that can build a PowerPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam ChildPreBasisFactory Types of child pre-basis factory
 * \param childPreBasisFactory Child pre-basis factory
 *
 * This overload will select the BasisFactory::BlockedInterleaved strategy.
 */
template<std::size_t k, class ChildPreBasisFactory>
auto power(ChildPreBasisFactory&& childPreBasisFactory)
{
  return Imp::PowerPreBasisFactory<k, BlockedInterleaved, ChildPreBasisFactory>(std::forward<ChildPreBasisFactory>(childPreBasisFactory));
}

} // end namespace BasisFactory

// Backward compatibility
namespace BasisBuilder {

  using namespace BasisFactory;

}


} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
