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
//   PowerNodeFactory
//   PowerNodeIndexSet
//
// The factory allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<class MI, class TP, class IMS, class SF, std::size_t C>
class PowerNodeIndexSet;



/**
 * \brief A factory for power bases
 *
 * This node factory represente a power of a given node factory.
 * Its node type is a PowerBasisNodes for the given subnode.
 *
 * \tparam MI  Type to be used for multi-indices
 * \tparam IMS An IndexMergingStrategy used to merge the global indices of the child factories
 * \tparam SF  The child factory
 * \tparam C   The exponent of the power node
 */
template<class MI, class IMS, class SF, std::size_t C>
class PowerNodeFactory
{
  static const std::size_t children = C;

  template<class, class, class, class, std::size_t>
  friend class PowerNodeIndexSet;

public:

  //! The child factory
  using SubFactory = SF;

  //! The grid view that the FE basis is defined on
  using GridView = typename SF::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Strategy used to merge the global indices of the child factories
  using IndexMergingStrategy = IMS;

  template<class TP>
  using SubNode = typename SubFactory::template Node<decltype(TypeTree::push_back(TP(), 0))>;

  template<class TP>
  using SubIndexSet = typename SubFactory::template IndexSet<decltype(TypeTree::push_back(TP(), 0))>;

  //! Template mapping root tree path to type of created tree node
  template<class TP>
  using Node = PowerBasisNode<size_type, TP, SubNode<TP>, children>;

  //! Template mapping root tree path to type of created tree node index set
  template<class TP>
  using IndexSet = PowerNodeIndexSet<MI, TP, IMS, SF, C>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, SubFactory::SizePrefix::max_size()+1>;

private:

  using SubMultiIndex = MI;

public:

  /**
   * \brief Constructor for given child factory objects
   *
   * The child factories will be stored as copies
   */
  template<class... SFArgs,
    disableCopyMove<PowerNodeFactory, SFArgs...> = 0,
    enableIfConstructible<SubFactory, SFArgs...> = 0>
  PowerNodeFactory(SFArgs&&... sfArgs) :
    subFactory_(std::forward<SFArgs>(sfArgs)...)
  {
    static_assert(models<Concept::NodeFactory<GridView>, SubFactory>(), "Subfactory passed to PowerNodeFactory does not model the NodeFactory concept.");
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    subFactory_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return subFactory_.gridView();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
    subFactory_.update(gv);
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
  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    auto node = Node<TP>(tp);
    for (std::size_t i=0; i<children; ++i)
      node.setChild(i, subFactory_.node(TypeTree::push_back(tp, i)));
    return node;
  }

  /**
   * \brief Create tree node index set with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * Create an index set suitable for the tree node obtained
   * by node(tp).
   */
  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
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

  size_type size(const SizePrefix& prefix, BasisBuilder::FlatInterleaved) const
  {
    // The root index size is the root index size of a single subnode
    // multiplied by the number of subnodes because we enumerate all
    // child indices in a row.
    if (prefix.size() == 0)
      return children*subFactory_.size({});

    // The first prefix entry refers to one of the (root index size)
    // subindex trees. Hence we have to first compute the corresponding
    // prefix entry for a single subnode subnode. The we can append
    // the other prefix entries unmodified, because the index tree
    // looks the same after the first level.
    typename SubFactory::SizePrefix subPrefix;
    subPrefix.push_back(prefix[0] / children);
    for(std::size_t i=1; i<prefix.size(); ++i)
      subPrefix.push_back(prefix[i]);
    return subFactory_.size(subPrefix);
  }

  size_type size(const SizePrefix& prefix, BasisBuilder::FlatLexicographic) const
  {
    // The size at the index tree root is the size of at the index tree
    // root of a single subnode multiplied by the number of subnodes
    // because we enumerate all child indices in a row.
    if (prefix.size() == 0)
      return children*subFactory_.size({});

    // The first prefix entry refers to one of the (root index size)
    // subindex trees. Hence we have to first compute the corresponding
    // prefix entry for a single subnode subnode. The we can append
    // the other prefix entries unmodified, because the index tree
    // looks the same after the first level.
    typename SubFactory::SizePrefix subPrefix;
    subPrefix.push_back(prefix[0] % children);
    for(std::size_t i=1; i<prefix.size(); ++i)
      subPrefix.push_back(prefix[i]);
    return subFactory_.size(subPrefix);
  }

  size_type size(const SizePrefix& prefix, BasisBuilder::BlockedLexicographic) const
  {
    if (prefix.size() == 0)
      return children;
    typename SubFactory::SizePrefix subPrefix;
    for(std::size_t i=1; i<prefix.size(); ++i)
      subPrefix.push_back(prefix[i]);
    return subFactory_.size(subPrefix);
  }

  size_type size(const SizePrefix& prefix, BasisBuilder::LeafBlockedInterleaved) const
  {
    if (prefix.size() == 0)
      return subFactory_.size();

    typename SubFactory::SizePrefix subPrefix;
    for(std::size_t i=0; i<prefix.size()-1; ++i)
      subPrefix.push_back(prefix[i]);

    size_type r = subFactory_.size(subPrefix);
    if (r==0)
      return 0;
    subPrefix.push_back(prefix.back());
    r = subFactory_.size(subPrefix);
    if (r==0)
      return children;
    return r;
  }

public:

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return subFactory_.dimension() * children;
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return subFactory_.maxNodeSize() * children;
  }

protected:
  SubFactory subFactory_;
};



template<class MI, class TP, class IMS, class SF, std::size_t C>
class PowerNodeIndexSet
{
  static const std::size_t children = C;

public:

  using SubFactory = SF;

  /** \brief The grid view that the FE space is defined on */
  using GridView = typename SF::GridView;
  using size_type = std::size_t;
  using IndexMergingStrategy = IMS;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PowerNodeFactory<MI, IMS, SF, C>;

  using Node = typename NodeFactory::template Node<TP>;

  using SubTreePath = typename TypeTree::Child<Node,0>::TreePath;

  using SubNodeIndexSet = typename NodeFactory::SubFactory::template IndexSet<SubTreePath>;

  PowerNodeIndexSet(const NodeFactory & nodeFactory) :
    nodeFactory_(&nodeFactory),
    subNodeIndexSet_(nodeFactory_->subFactory_.template indexSet<SubTreePath>())
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

  template<typename It>
  It indices(It multiIndices, BasisBuilder::FlatInterleaved) const
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
  It indices(It multiIndices, BasisBuilder::FlatLexicographic) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    size_type firstIndexEntrySize = nodeFactory_->subFactory_.size({});
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

  static const void multiIndexPushFront(MultiIndex& M, size_type M0)
  {
    M.resize(M.size()+1);
    for(std::size_t i=M.size()-1; i>0; --i)
      M[i] = M[i-1];
    M[0] = M0;
  }

  template<typename It>
  It indices(It multiIndices, BasisBuilder::BlockedLexicographic) const
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
  It indices(It multiIndices, BasisBuilder::LeafBlockedInterleaved) const
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

private:
  const NodeFactory* nodeFactory_;
  SubNodeIndexSet subNodeIndexSet_;
  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t k, class IndexMergingStrategy, class SubFactoryTag>
struct PowerNodeFactoryBuilder
{
  static const bool isBlocked = std::is_same<IndexMergingStrategy,BlockedLexicographic>::value or std::is_same<IndexMergingStrategy,LeafBlockedInterleaved>::value;

  static const std::size_t requiredMultiIndexSize=SubFactoryTag::requiredMultiIndexSize + (std::size_t)(isBlocked);

  template<class MultiIndex, class GridView>
  auto build(const GridView& gridView)
    -> PowerNodeFactory<MultiIndex,  IndexMergingStrategy, decltype(SubFactoryTag().template build<MultiIndex, GridView>(std::declval<GridView>())), k>
  {
    return {SubFactoryTag().template build<MultiIndex, GridView>(gridView)};
  }
};

} // end namespace BasisBuilder::Imp



/**
 * \brief Create a factory builder that can build a PowerNodeFactory
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam SubFactoryTag Types of child factory builder and IndexMergingStrategy type
 * \tparam IndexMergingStrategy An IndexMergingStrategy type
 * \param tag Child factory builder objects and an IndexMergingStrategy
 * \param ims IndexMergingStrategy to be used
 *
 * This overload can be used to explicitly supply an IndexMergingStrategy.
 */
template<std::size_t k, class SubFactoryTag, class IndexMergingStrategy>
Imp::PowerNodeFactoryBuilder<k, IndexMergingStrategy, SubFactoryTag>
  power(SubFactoryTag&& tag, const IndexMergingStrategy& ims)
{
  return{};
}

/**
 * \brief Create a factory builder that can build a PowerNodeFactory
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam SubFactoryTag Types of child factory builder and IndexMergingStrategy type
 * \param tag Child factory builder objects and an IndexMergingStrategy
 *
 * This overload will select the BasisBuilder::LeafBlockedInterleaved strategy.
 */
template<std::size_t k, class SubFactoryTag>
Imp::PowerNodeFactoryBuilder<k, LeafBlockedInterleaved, SubFactoryTag>
  power(SubFactoryTag&& tag)
{
  return{};
}

} // end namespace BasisBuilder



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
