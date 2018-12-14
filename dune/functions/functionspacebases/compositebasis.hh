// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH

#include <tuple>
#include <utility>

#include <dune/common/std/utility.hh>
#include <dune/common/std/apply.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/tuplevector.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/utility.hh>

#include <dune/functions/common/staticforloop.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the composite bases. It contains
//
//   CompositePreBasis
//   CompositeNodeIndexSet
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************


template<class MI, class TP, class IT, class... SPB>
class CompositeNodeIndexSet;

/**
 * \brief A pre-basis for composite bases
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * This pre-basis represente a composition of several given pre-bases.
 * Its node type is a CompositeBasisNodes for the given subnodes.
 *
 * \tparam MI  Type to be used for global multi-indices
 * \tparam IMS An IndexMergingStrategy used to merge the global indices of the child pre-bases
 * \tparam SPB  The child pre-bases
 */
template<class MI, class IMS, class... SPB>
class CompositePreBasis
{
public:

  //! Tuple of child pre-bases
  using SubPreBases = std::tuple<SPB...>;

  //! The grid view that the FE basis is defined on
  using GridView = typename std::tuple_element_t<0, SubPreBases>::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Strategy used to merge the global indices of the child pre-bases
  using IndexMergingStrategy = IMS;

protected:
  static const std::size_t children = sizeof...(SPB);

  template<class, class, class, class...>
  friend class CompositeNodeIndexSet;

  using ChildIndices = std::make_index_sequence<children>;

  template<class TP, class Indices>
  struct Types;

  template<class TP, std::size_t... indices>
  struct Types<TP, Dune::Std::index_sequence<indices...>>
  {
    template<std::size_t i>
    using SubTreePath = decltype(TypeTree::push_back(TP(), Dune::index_constant<i>()));

    template<std::size_t i>
    using SubNode = typename std::tuple_element_t<i, SubPreBases>::template Node<SubTreePath<i>>;

    template<std::size_t i>
    using SubIndexSet = typename std::tuple_element_t<i, SubPreBases>::template IndexSet<SubTreePath<i>>;

    using SubIndexSets = std::tuple<SubIndexSet<indices>...>;

    using Node = CompositeBasisNode<SubNode<indices>...>;
  };

  template<class TP>
  using SubIndexSets = typename Types<TP, ChildIndices>::SubIndexSets;

public:

  //! Template mapping index of child to its pre-basis type
  template<std::size_t k>
  using SubPreBasis = typename std::tuple_element<k, std::tuple<SPB...>>::type;

  //! Template mapping root tree path to type of created tree node
  template<class TP>
  using Node = typename Types<TP, ChildIndices>::Node;

  //! Template mapping root tree path to type of created tree node index set
  template<class TP>
  using IndexSet = CompositeNodeIndexSet<MI, TP, IMS, SPB...>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, MultiIndex::max_size()>;

  /**
   * \brief Constructor for given child pre-basis objects
   *
   * The child pre-basis will be stored as copies
   */
  template<class... SFArgs,
    disableCopyMove<CompositePreBasis, SFArgs...> = 0,
    enableIfConstructible<std::tuple<SPB...>, SFArgs...> = 0>
  CompositePreBasis(SFArgs&&... sfArgs) :
    subPreBases_(std::forward<SFArgs>(sfArgs)...)
  {
    Hybrid::forEach(subPreBases_, [&](const auto& subPreBasis){
      static_assert(models<Concept::PreBasis<GridView>, std::decay_t<decltype(subPreBasis)>>(), "Subprebases passed to CompositePreBasis does not model the PreBasis concept.");
    });
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      Hybrid::elementAt(subPreBases_, i).initializeIndices();
    });
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return std::get<0>(subPreBases_).gridView();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      Hybrid::elementAt(subPreBases_, i).update(gv);
    });
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
    auto node = Node<TP>();
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      node.setChild(Hybrid::elementAt(subPreBases_, i).node(TypeTree::push_back(tp, i)), i);
    });
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
    return IndexSet<TP>{*this,
      unpackIntegerSequence([&](auto... i) {
        return std::make_tuple(Hybrid::elementAt(subPreBases_, i).template indexSet<decltype(TypeTree::push_back(TP(), i))>()...);
      }, ChildIndices())};
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

  size_type size(const SizePrefix& prefix, BasisFactory::BlockedLexicographic) const
  {
    if (prefix.size() == 0)
      return children;

    return Hybrid::switchCases(ChildIndices(), prefix[0], [&] (auto i) {
      const auto& subPreBasis = std::get<i.value>(subPreBases_);
      typename std::decay<decltype(subPreBasis)>::type::SizePrefix subPrefix;
      for(std::size_t i=1; i<prefix.size(); ++i)
        subPrefix.push_back(prefix[i]);
      return subPreBasis.size(subPrefix);
    }, []() {
      return size_type(0);
    });
  }

  struct Lambda_size_flat_sizeInSubtree
  {
    template<class I, class PB>
    size_type operator()(const I&, const PB& subPreBases, const SizePrefix& prefix, size_type& shiftedFirst, size_type& r)
    {
      using SubPreBasis = typename std::tuple_element<I::value, PB>::type;
      const SubPreBasis& subPreBasis = std::get<I::value>(subPreBases);
      if (shiftedFirst < subPreBasis.size())
      {
        typename SubPreBasis::SizePrefix subPrefix;
        subPrefix.push_back(shiftedFirst);
        for(std::size_t i=1; i<prefix.size(); ++i)
          subPrefix.push_back(prefix[i]);
        r = subPreBasis.size(subPrefix);
        return true;
      }
      shiftedFirst -= subPreBasis.size();
      return false;
    }
  };

  size_type size(const SizePrefix& prefix, BasisFactory::FlatLexicographic) const
  {
    size_type r = 0;
    if (prefix.size() == 0)
      Hybrid::forEach(ChildIndices(), [&](auto i) {
        r += Hybrid::elementAt(subPreBases_, i).size();
      });
    else {
      size_type shiftedFirst = prefix[0];
      staticFindInRange<0, sizeof...(SPB)>(Lambda_size_flat_sizeInSubtree(), subPreBases_, prefix, shiftedFirst, r);
    }
    return r;
  }

public:

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    size_type r=0;
    // Accumulate dimension() for all subprebases
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      r += Hybrid::elementAt(subPreBases_, i).dimension();
    });
    return r;
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    size_type r=0;
    // Accumulate maxNodeSize() for all subprebases
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      r += Hybrid::elementAt(subPreBases_, i).maxNodeSize();
    });
    return r;
  }

protected:
  std::tuple<SPB...> subPreBases_;
};



template<class MI, class TP, class IMS, class... SPB>
class CompositeNodeIndexSet
{
  static const std::size_t children = sizeof...(SPB);

  using ChildIndices = std::make_index_sequence<children>;
public:

  using PreBasis = CompositePreBasis<MI, IMS, SPB...>;

  using GridView = typename PreBasis::GridView;
  using size_type = std::size_t;
  using IndexMergingStrategy = IMS;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using Node = typename PreBasis::template Node<TP>;

  using SubIndexSets = typename PreBasis::template SubIndexSets<TP>;

  CompositeNodeIndexSet(const PreBasis & preBasis, SubIndexSets&& subNodeIndexSets_) :
    preBasis_(&preBasis),
    subNodeIndexSetTuple_(std::move(subNodeIndexSets_)),
    node_(nullptr)
  {}

  void bind(const Node& node)
  {
    node_ = &node;
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      Hybrid::elementAt(subNodeIndexSetTuple_, i).bind(node.child(i));
    });
  }

  void unbind()
  {
    node_ = nullptr;
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      Hybrid::elementAt(subNodeIndexSetTuple_, i).unbind();
    });
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
  It indices(It multiIndices, BasisFactory::FlatLexicographic) const
  {
    size_type firstComponentOffset = 0;
    // Loop over all children
    Hybrid::forEach(ChildIndices(), [&](auto child){
      const auto& subNodeIndexSet = Hybrid::elementAt(subNodeIndexSetTuple_, child);
      const auto& subPreBasis = Hybrid::elementAt(preBasis_->subPreBases_, child);
      size_type subTreeSize = subNodeIndexSet.size();
      // Fill indices for current child into index buffer starting from current
      // buffer position and shift first index component of any index for current
      // child by suitable offset to get lexicographic indices.
      subNodeIndexSet.indices(multiIndices);
      for (std::size_t i = 0; i<subTreeSize; ++i)
        multiIndices[i][0] += firstComponentOffset;
      // Increment offset by the size for first index component of the current child
      firstComponentOffset += subPreBasis.size({});
      // Increment buffer iterator by the number of indices processed for current child
      multiIndices += subTreeSize;
    });
    return multiIndices;
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
    // Loop over all children
    Hybrid::forEach(ChildIndices(), [&](auto child){
      const auto& subNodeIndexSet = Hybrid::elementAt(subNodeIndexSetTuple_, child);
      size_type subTreeSize = subNodeIndexSet.size();
      // Fill indices for current child into index buffer starting from current position
      subNodeIndexSet.indices(multiIndices);
      // Insert child index before first component of all indices of current child.
      for (std::size_t i = 0; i<subTreeSize; ++i)
        this->multiIndexPushFront(multiIndices[i], child);
      // Increment buffer iterator by the number of indices processed for current child
      multiIndices += subTreeSize;
    });
    return multiIndices;
  }


private:
  const PreBasis* preBasis_;
  SubIndexSets subNodeIndexSetTuple_;
  const Node* node_;
};



namespace BasisFactory {

namespace Imp {

template<class ST0>
constexpr std::size_t maxHelper(ST0&& i0)
{
  return i0;
}

template<class ST0, class... ST>
constexpr std::size_t maxHelper(ST0&& i0, ST&&... i)
{
  return (i0 > maxHelper(i...)) ? i0 : maxHelper(i...);
}

template<class IndexMergingStrategy, class... ChildPreBasisFactory>
class CompositePreBasisFactory
{
  static const bool isBlocked = std::is_same<IndexMergingStrategy,BlockedLexicographic>::value or std::is_same<IndexMergingStrategy,BlockedInterleaved>::value;

  static const std::size_t maxChildIndexSize = maxHelper(ChildPreBasisFactory::requiredMultiIndexSize...);

  template<class MultiIndex, class GridView, class... ChildPreBasis>
  auto makePreBasisFromChildPreBases(const GridView&, ChildPreBasis&&... childPreBasis) const
  {
    return CompositePreBasis<MultiIndex, IndexMergingStrategy, std::decay_t<ChildPreBasis>...>(std::forward<ChildPreBasis>(childPreBasis)...);
  }

public:

  static const std::size_t requiredMultiIndexSize = isBlocked ? (maxChildIndexSize+1) : maxChildIndexSize;

  CompositePreBasisFactory(const ChildPreBasisFactory&... childPreBasisFactory) :
    childPreBasisFactories_(childPreBasisFactory...)
  {}

  CompositePreBasisFactory(ChildPreBasisFactory&&... childPreBasisFactory) :
    childPreBasisFactories_(std::move(childPreBasisFactory)...)
  {}

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    // Use Std::apply to unpack the tuple childPreBasisFactories_
    return Std::apply([&](const auto&... childPreBasisFactory) {
        return this->makePreBasisFromChildPreBases<MultiIndex>(gridView, childPreBasisFactory.template makePreBasis<MultiIndex>(gridView)...);
      }, childPreBasisFactories_);
  }

private:
  std::tuple<ChildPreBasisFactory...> childPreBasisFactories_;
};

} // end namespace BasisFactory::Imp



/**
 * \brief Create a factory builder that can build a CompositePreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam Args Types of child factory builders and IndexMergingStrategy type
 * \param args Child factory builder objects and an IndexMergingStrategy
 *
 * This is the overload used if the last argument is an IndexMergingStrategy.
 */
template<
  typename... Args,
  std::enable_if_t<Concept::isIndexMergingStrategy<typename LastType<Args...>::type>(),int> = 0>
auto composite(Args&&... args)
{
  // We have to separate the last entry which is the IndexMergingStrategy
  // and the preceding ones, which are the ChildPreBasisFactories

  using ArgTuple = std::tuple<std::decay_t<Args>...>;

  // Compute number of children and index of the IndexMergingStrategy argument
  constexpr std::size_t children = Dune::SizeOf<Args...>::value-1;

  // Use last type as IndexMergingStrategy
  using IndexMergingStrategy = std::tuple_element_t<children, ArgTuple>;

  // Index sequence for all but the last entry for partial tuple unpacking
  auto childIndices = std::make_index_sequence<children>{};

  // Unpack tuple only for those entries related to children
  return applyPartial([](auto&&... childPreBasisFactory){
    return Imp::CompositePreBasisFactory<IndexMergingStrategy, std::decay_t<decltype(childPreBasisFactory)>...>(std::forward<decltype(childPreBasisFactory)>(childPreBasisFactory)...);
  },
  std::forward_as_tuple(std::forward<Args>(args)...),
  childIndices);
}

/**
 * \brief Create a factory builder that can build a CompositePreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam Args Types of child factory builders
 * \param args Child factory builder objects
 *
 * This is the overload used if no IndexMergingStrategy is supplied.
 * In this case the BasisFactory::BlockedLexicographic strategy is used.
 */
template<
  typename... Args,
  std::enable_if_t<not Concept::isIndexMergingStrategy<typename LastType<Args...>::type>(),int> = 0>
auto composite(Args&&... args)
{
  return Imp::CompositePreBasisFactory<BasisFactory::BlockedLexicographic, std::decay_t<Args>...>(std::forward<Args>(args)...);
}

} // end namespace BasisFactory

// Backward compatibility
namespace BasisBuilder {

  using namespace BasisFactory;

}



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH
