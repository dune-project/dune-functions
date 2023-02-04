// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH

#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/tuplevector.hh>

#include <dune/functions/common/staticforloop.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the composite bases. It contains
//
//   CompositePreBasis
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************


/**
 * \brief A pre-basis for composite bases
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * This pre-basis represents a composition of several given pre-bases.
 * Its node type is a CompositeBasisNodes for the given subnodes.
 *
 * \tparam IMS An IndexMergingStrategy used to merge the global indices of the child pre-bases
 * \tparam SPB  The child pre-bases
 */
template<class IMS, class... SPB>
class CompositePreBasis
{
  static const bool isBlocked = std::is_same_v<IMS,BasisFactory::BlockedLexicographic> or std::is_same_v<IMS,BasisFactory::BlockedInterleaved>;
public:

  //! Tuple of child pre-bases
  using SubPreBases = std::tuple<SPB...>;

  //! Export individual child pre-bases by index
  template<std::size_t i>
  using SubPreBasis = std::tuple_element_t<i, SubPreBases>;

  //! The grid view that the FE basis is defined on
  using GridView = typename std::tuple_element_t<0, SubPreBases>::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Strategy used to merge the global indices of the child pre-bases
  using IndexMergingStrategy = IMS;

protected:
  static const std::size_t children = sizeof...(SPB);

  using ChildIndices = std::make_index_sequence<children>;

public:

  //! Template mapping root tree path to type of created tree node
  using Node = CompositeBasisNode<typename SPB::Node...>;

  static constexpr size_type maxMultiIndexSize = std::max({SPB::maxMultiIndexSize...}) + isBlocked;
  static constexpr size_type minMultiIndexSize = std::min({SPB::minMultiIndexSize...}) + isBlocked;
  static constexpr size_type multiIndexBufferSize = std::max({SPB::multiIndexBufferSize...}) + isBlocked;

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

  /**
   * \brief Constructor for given GridView
   *
   * This constructor is only available if all child pre-bases are constructible
   * from the grid view.
   */
  template<class GV,
    std::enable_if_t<std::conjunction_v<
      std::bool_constant<(children > 1)>,    // Avoid ambiguous constructor if there's only one child
      std::is_same<GV, GridView>,
      std::is_constructible<SPB, GridView>...
    >, int> = 0>
  CompositePreBasis(const GV& gv) :
    subPreBases_(SPB(gv)...)
  {
    Hybrid::forEach(subPreBases_, [&](const auto& subPreBasis){
      static_assert(models<Concept::PreBasis<GridView>, std::decay_t<decltype(subPreBasis)>>(), "Subprebases passed to CompositePreBasis does not model the PreBasis concept.");
    });
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      this->subPreBasis(i).initializeIndices();
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
      this->subPreBasis(i).update(gv);
    });
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    auto node = Node{};
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      node.setChild(this->subPreBasis(i).makeNode(), i);
    });
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
  size_type size(const SizePrefix& prefix, BasisFactory::BlockedLexicographic) const
  {
    if (prefix.size() == 0)
      return children;

    return Hybrid::switchCases(ChildIndices(), prefix[0], [&] (auto i) {
      SizePrefix subPrefix;
      for(std::size_t i=1; i<prefix.size(); ++i)
        subPrefix.push_back(prefix[i]);
      return this->subPreBasis(i).size(subPrefix);
    }, []() {
      return size_type(0);
    });
  }

  template<class SizePrefix>
  size_type size(const SizePrefix& prefix, BasisFactory::FlatLexicographic) const
  {
    size_type result = 0;
    if (prefix.size() == 0)
      Hybrid::forEach(ChildIndices(), [&](auto i) {
        result += this->subPreBasis(i).size();
      });
    else {
      size_type shiftedFirstDigit = prefix[0];
      staticFindInRange<0, children>([&](auto i) {
          auto firstDigitSize = this->subPreBasis(i).size();
          if (shiftedFirstDigit < firstDigitSize)
          {
            SizePrefix subPrefix;
            subPrefix.push_back(shiftedFirstDigit);
            for(std::size_t i=1; i<prefix.size(); ++i)
              subPrefix.push_back(prefix[i]);
            result = this->subPreBasis(i).size(subPrefix);
            return true;
          }
          shiftedFirstDigit -= firstDigitSize;
          return false;
        });
    }
    return result;
  }

public:

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    size_type r=0;
    // Accumulate dimension() for all subprebases
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      r += this->subPreBasis(i).dimension();
    });
    return r;
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    size_type r=0;
    // Accumulate maxNodeSize() for all subprebases
    Hybrid::forEach(ChildIndices(), [&](auto i) {
      r += this->subPreBasis(i).maxNodeSize();
    });
    return r;
  }

  //! Const access to the stored prebasis of the factor in the power space
  template<std::size_t i>
  const SubPreBasis<i>& subPreBasis(Dune::index_constant<i> = {}) const
  {
    return std::get<i>(subPreBases_);
  }

  //! Mutable access to the stored prebasis of the factor in the power space
  template<std::size_t i>
  SubPreBasis<i>& subPreBasis(Dune::index_constant<i> = {})
  {
    return std::get<i>(subPreBases_);
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it) const
  {
    return indices(node, it, IndexMergingStrategy{});
  }

private:

  template<typename It>
  It indices(const Node& node, It multiIndices, BasisFactory::FlatLexicographic) const
  {
    size_type firstComponentOffset = 0;
    // Loop over all children
    Hybrid::forEach(ChildIndices(), [&](auto child){
      size_type subTreeSize = node.child(child).size();
      // Fill indices for current child into index buffer starting from current
      // buffer position and shift first index component of any index for current
      // child by suitable offset to get lexicographic indices.
      subPreBasis(child).indices(node.child(child), multiIndices);
      for (std::size_t i = 0; i<subTreeSize; ++i)
        multiIndices[i][0] += firstComponentOffset;
      // Increment offset by the size for first index component of the current child
      firstComponentOffset += subPreBasis(child).size();
      // Increment buffer iterator by the number of indices processed for current child
      multiIndices += subTreeSize;
    });
    return multiIndices;
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
    // Loop over all children
    Hybrid::forEach(ChildIndices(), [&](auto child){
      size_type subTreeSize = node.child(child).size();
      // Fill indices for current child into index buffer starting from current position
      subPreBasis(child).indices(node.child(child), multiIndices);
      // Insert child index before first component of all indices of current child.
      for (std::size_t i = 0; i<subTreeSize; ++i)
        this->multiIndexPushFront(multiIndices[i], child);
      // Increment buffer iterator by the number of indices processed for current child
      multiIndices += subTreeSize;
    });
    return multiIndices;
  }

  std::tuple<SPB...> subPreBases_;
};



namespace BasisFactory {

namespace Imp {

template<class IndexMergingStrategy, class... ChildPreBasisFactory>
class CompositePreBasisFactory
{

  template<class GridView, class... ChildPreBasis>
  auto makePreBasisFromChildPreBases(const GridView&, ChildPreBasis&&... childPreBasis) const
  {
    return CompositePreBasis<IndexMergingStrategy, std::decay_t<ChildPreBasis>...>(std::forward<ChildPreBasis>(childPreBasis)...);
  }

public:

  CompositePreBasisFactory(const ChildPreBasisFactory&... childPreBasisFactory) :
    childPreBasisFactories_(childPreBasisFactory...)
  {}

  CompositePreBasisFactory(ChildPreBasisFactory&&... childPreBasisFactory) :
    childPreBasisFactories_(std::move(childPreBasisFactory)...)
  {}

  template<class GridView>
  auto operator()(const GridView& gridView) const
  {
    // Use std::apply to unpack the tuple childPreBasisFactories_
    return std::apply([&](const auto&... childPreBasisFactory) {
        return this->makePreBasisFromChildPreBases(gridView, childPreBasisFactory(gridView)...);
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
