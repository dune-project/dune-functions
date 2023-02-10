// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>



namespace Dune {
namespace Functions {


// *****************************************************************************
// This is the reusable part of the power bases. It contains
//
//   PowerPreBasis
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

/**
 * \brief A pre-basis for power bases
 *
 * This pre-basis represents a power of a given pre-basis.
 * Its node type is a PowerBasisNodes for the given subnode.
 *
 * \tparam IMS An IndexMergingStrategy used to merge the global indices of the child factories
 * \tparam SPB  The child pre-basis
 * \tparam C   The exponent of the power node
 */
template<class IMS, class SPB, std::size_t C>
class PowerPreBasis
{
  static const std::size_t children = C;
  static const bool isBlocked = std::is_same_v<IMS,BasisFactory::BlockedLexicographic> or std::is_same_v<IMS,BasisFactory::BlockedInterleaved>;

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

  //! Template mapping root tree path to type of created tree node
  using Node = PowerBasisNode<SubNode, children>;

  static constexpr size_type maxMultiIndexSize = SubPreBasis::maxMultiIndexSize + isBlocked;
  static constexpr size_type minMultiIndexSize = SubPreBasis::minMultiIndexSize + isBlocked;
  static constexpr size_type multiIndexBufferSize = SubPreBasis::multiIndexBufferSize + isBlocked;

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

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return size(Dune::ReservedVector<size_type, multiIndexBufferSize>{});
  }

  //! Return number of possible values for next position in multi index
  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
  {
    return size(prefix, IndexMergingStrategy{});
  }

private:

  template<class SizePrefix>
  size_type size(SizePrefix prefix, BasisFactory::FlatInterleaved) const
  {
    // The root index size is the root index size of a single subnode
    // multiplied by the number of subnodes, because we enumerate all
    // child indices in a row.
    if (prefix.size() == 0)
      return children*subPreBasis_.size();

    // The FlatInterleaved index merging strategy only changes the first
    // index digit. Hence, we have to reconstruct the corresponding digit
    // for the subtree and can then return the corresponding size of the subtree.
    prefix[0] = prefix[0] / children;
    return subPreBasis_.size(prefix);
  }

  template<class SizePrefix>
  size_type size(SizePrefix prefix, BasisFactory::FlatLexicographic) const
  {
    // The size at the index tree root is the size of at the index tree
    // root of a single subnode multiplied by the number of subnodes,
    // because we enumerate all child indices in a row.
    if (prefix.size() == 0)
      return children*subPreBasis_.size();

    // The first prefix entry refers to one of the (root index size)
    // subindex trees. Hence, we have to first compute the corresponding
    // prefix entry for a single subnode subnode. Then we can append
    // the other prefix entries unmodified, because the index tree
    // looks the same after the first level.

    // The FlatLexicographic index merging strategy only changes the first
    // index digit. Hence, we have to reconstruct the corresponding digit
    // for the subtree and can then return the corresponding size of the subtree.
    prefix[0] = prefix[0] % subPreBasis_.size();
    return subPreBasis_.size(prefix);
  }

  template<class MultiIndex>
  static void multiIndexPopFront(MultiIndex& M)
  {
    for(std::size_t i=0; i<M.size()-1; ++i)
      M[i] = M[i+1];
    M.resize(M.size()-1);
  }

  template<class SizePrefix>
  size_type size(SizePrefix prefix, BasisFactory::BlockedLexicographic) const
  {
    if (prefix.size() == 0)
      return children;
    multiIndexPopFront(prefix);
    return subPreBasis_.size(prefix);
  }

  template<class SizePrefix>
  size_type size(SizePrefix prefix, BasisFactory::BlockedInterleaved) const
  {
    if (prefix.size() == 0)
      return subPreBasis_.size();

    // Remember last index, remove it and check if the remaining
    // prefix refers to a leaf in the subPreBasis index tree.
    // If yes, then the full prefix must also refer to a
    // leaf in the merged index tree. If not, then restore the full
    // prefix and proceed.
    auto tail = prefix.back();
    prefix.pop_back();
    if (subPreBasis_.size(prefix) == 0)
      return 0;
    prefix.push_back(tail);

    // Now check if the full prefix refers to a leaf in the subPreBasis
    // index tree.
    // If yes, then it has exactly 'children' appended children in the subtree.
    // If not, then the index tree looks the same in the merged subtree and we
    // can forward the result.
    auto subSize = subPreBasis_.size(prefix);
    if (subSize == 0)
      return children;
    return subSize;
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

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it) const
  {
    return indices(node, it, IndexMergingStrategy{});
  }

private:

  template<typename It>
  It indices(const Node& node, It multiIndices, BasisFactory::FlatInterleaved) const
  {
    using namespace Dune::Indices;
    size_type subTreeSize = node.child(_0).size();
    // Fill indices for first child at the beginning.
    auto next = subPreBasis().indices(node.child(_0), multiIndices);
    // Multiply first component of all indices for first child by
    // number of children to stretch the index range for interleaving.
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
  It indices(const Node& node, It multiIndices, BasisFactory::FlatLexicographic) const
  {
    using namespace Dune::Indices;
    size_type subTreeSize = node.child(_0).size();
    size_type firstIndexEntrySize = subPreBasis().size();
    // Fill indices for first child at the beginning.
    auto next = subPreBasis().indices(node.child(_0), multiIndices);
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

  template<class MultiIndex>
  static void multiIndexPushFront(MultiIndex& M, size_type M0)
  {
    M.resize(M.size()+1);
    for(std::size_t i=M.size()-1; i>0; --i)
      M[i] = M[i-1];
    M[0] = M0;
  }

  template<typename It>
  It indices(const Node& node, It multiIndices, BasisFactory::BlockedLexicographic) const
  {
    using namespace Dune::Indices;
    size_type subTreeSize = node.child(_0).size();
    // Fill indices for first child at the beginning.
    auto next = subPreBasis().indices(node.child(_0), multiIndices);
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
  It indices(const Node& node, It multiIndices, BasisFactory::BlockedInterleaved) const
  {
    using namespace Dune::Indices;
    size_type subTreeSize = node.child(_0).size();
    // Fill indices for first child at the beginning.
    auto next = subPreBasis().indices(node.child(_0), multiIndices);
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

  SubPreBasis subPreBasis_;
};



namespace BasisFactory {

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
auto power(ChildPreBasisFactory&& childPreBasisFactory, const IndexMergingStrategy&)
{
  return [childPreBasisFactory](const auto& gridView) {
    auto childPreBasis = childPreBasisFactory(gridView);
    return PowerPreBasis<IndexMergingStrategy, decltype(childPreBasis), k>(std::move(childPreBasis));
  };
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
  return [childPreBasisFactory](const auto& gridView) {
    auto childPreBasis = childPreBasisFactory(gridView);
    return PowerPreBasis<BlockedInterleaved, decltype(childPreBasis), k>(std::move(childPreBasis));
  };
}

} // end namespace BasisFactory

// Backward compatibility
namespace BasisBuilder {

  using namespace BasisFactory;

}


} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
