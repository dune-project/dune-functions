// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH

#include <dune/common/reservedvector.hh>

#include <dune/typetree/powernode.hh>
#include <dune/typetree/utility.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/staticforloop.hh>
#include <dune/functions/functionspacebases/basistags.hh>



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

template<class MI, class TP, class IT, class SF, std::size_t C>
class PowerNodeIndexSet;



/**
 * \brief A factory for power bases
 *
 * This node factory represente a power of a given node factory.
 * Its node type is a PowerBasisNodes for the given subnode.
 *
 * \tparam MI Type to be used for multi-indices
 * \tparam IT A tag describing how global indices are build
 * \tparam SF The subnode factory
 * \tparam C The exponent of the power node
 */
template<class MI, class IT, class SF, std::size_t C>
class PowerNodeFactory
{
  static const std::size_t children = C;

  template<class, class, class, class, std::size_t>
  friend class PowerNodeIndexSet;

public:

  using SubFactory = SF;

  /** \brief The grid view that the FE space is defined on */
  using GridView = typename SF::GridView;
  using size_type = typename SF::size_type;
  using IndexTag = IT;

  template<class TP>
  using SubNode = typename SubFactory::template Node<decltype(TypeTree::push_back(TP(), 0))>;

  template<class TP>
  using SubIndexSet = typename SubFactory::template IndexSet<decltype(TypeTree::push_back(TP(), 0))>;

  template<class TP>
  using Node = PowerBasisNode<size_type, TP, SubNode<TP>, children>;

  template<class TP>
  using IndexSet = PowerNodeIndexSet<MI, TP, IT, SF, C>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, SubFactory::SizePrefix::max_size()+1>;

private:

  using SubMultiIndex = MI;

public:

  /** \brief Constructor for a given grid view object */

  template<class... SFArgs>
  PowerNodeFactory(SFArgs&&... sfArgs) :
    subFactory_(std::forward<SFArgs>(sfArgs)...)
  {}


  void initializeIndices()
  {
    subFactory_.initializeIndices();
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return subFactory_.gridView();
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    auto node = Node<TP>(tp);
    for(int i=0; i<children; ++i)
      node.setChild(i, std::make_shared<SubNode<TP> >(subFactory_.node(TypeTree::push_back(tp, i))));
    return node;
  }

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

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return size(prefix, IndexTag());
  }

  size_type size(const SizePrefix& prefix, BasisTags::InterleafedIndex) const
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

  size_type size(const SizePrefix& prefix, BasisTags::FlatIndex) const
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

  size_type size(const SizePrefix& prefix, BasisTags::BlockedIndex) const
  {
    if (prefix.size() == 0)
      return children;
    typename SubFactory::SizePrefix subPrefix;
    for(std::size_t i=1; i<prefix.size(); ++i)
      subPrefix.push_back(prefix[i]);
    return subFactory_.size(subPrefix);
  }

  size_type size(const SizePrefix& prefix, BasisTags::LeafBlockedIndex) const
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

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return subFactory_.dimension() * children;
  }

  size_type maxNodeSize() const
  {
    return subFactory_.maxNodeSize() * children;
  }

protected:
  SubFactory subFactory_;
};



template<class MI, class TP, class IT, class SF, std::size_t C>
class PowerNodeIndexSet
{
  static const std::size_t children = C;

public:

  using SubFactory = SF;

  /** \brief The grid view that the FE space is defined on */
  using GridView = typename SF::GridView;
  using size_type = typename SF::size_type;
  using IndexTag = IT;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PowerNodeFactory<MI, IT, SF, C>;

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

  MultiIndex index(const size_type& localIndex) const
  {
    return index(localIndex, IndexTag());
  }


  MultiIndex index(const size_type& localIndex, BasisTags::InterleafedIndex) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    size_type subLocalIndex = localIndex % subTreeSize;
    size_type component = localIndex / subTreeSize;

    MultiIndex mi = subNodeIndexSet_.index(subLocalIndex);
    mi[0] = mi[0]*children+component;

    return mi;
  }

  MultiIndex index(const size_type& localIndex, BasisTags::FlatIndex) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    size_type subLocalIndex = localIndex % subTreeSize;
    size_type component = localIndex / subTreeSize;

    size_type firstLevelSize = nodeFactory_->subFactory_.size({});

    MultiIndex mi = subNodeIndexSet_.index(subLocalIndex);
    mi[0] += component*firstLevelSize;

    return mi;
  }

  MultiIndex index(const size_type& localIndex, BasisTags::BlockedIndex) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    size_type subLocalIndex = localIndex % subTreeSize;
    size_type component = localIndex / subTreeSize;

    auto subTreeMi = subNodeIndexSet_.index(subLocalIndex);

    MultiIndex mi;
    mi.push_back(component);
    for(std::size_t i=0; i<subTreeMi.size(); ++i)
      mi.push_back(subTreeMi[i]);
    return mi;
  }

  MultiIndex index(const size_type& localIndex, BasisTags::LeafBlockedIndex) const
  {
    using namespace Dune::TypeTree::Indices;
    size_type subTreeSize = node_->child(_0).size();
    size_type subLocalIndex = localIndex % subTreeSize;
    size_type component = localIndex / subTreeSize;

    auto subTreeMi = subNodeIndexSet_.index(subLocalIndex);

    MultiIndex mi;
    for(std::size_t i=0; i<subTreeMi.size(); ++i)
      mi.push_back(subTreeMi[i]);
    mi.push_back(component);
    return mi;
  }

private:
  const NodeFactory* nodeFactory_;
  SubNodeIndexSet subNodeIndexSet_;
  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t k, class IndexTag, class SubFactoryTag>
struct PowerNodeFactoryBuilder
{
  static const bool isBlocked = std::is_same<IndexTag,BasisTags::BlockedIndex>::value or std::is_same<IndexTag,BasisTags::LeafBlockedIndex>::value;

  static const std::size_t requiredMultiIndexSize=SubFactoryTag::requiredMultiIndexSize + (std::size_t)(isBlocked);

  template<class MultiIndex, class GridView, class size_type=std::size_t>
  auto build(const GridView& gridView)
    -> PowerNodeFactory<MultiIndex,  IndexTag, decltype(SubFactoryTag().template build<MultiIndex, GridView, size_type>(std::declval<GridView>())), k>
  {
    return {SubFactoryTag().build<MultiIndex, GridView, size_type>(gridView)};
  }
};

} // end namespace BasisBuilder::Imp

template<std::size_t k, class SubFactoryTag, class IndexTag>
Imp::PowerNodeFactoryBuilder<k, IndexTag, SubFactoryTag> power(SubFactoryTag&& tag, const IndexTag&)
{
  return{};
}

template<std::size_t k, class SubFactoryTag>
Imp::PowerNodeFactoryBuilder<k, BasisTags::LeafBlockedIndex, SubFactoryTag> power(SubFactoryTag&& tag)
{
  return{};
}

} // end namespace BasisBuilder



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
