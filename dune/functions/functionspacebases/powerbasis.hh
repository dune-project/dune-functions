// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/dynamicpowerbasis.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>
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
 * \tparam C    The exponent of the power node
 */
template<class IMS, class SPB, std::size_t C>
class PowerPreBasis :
    public DynamicPowerPreBasis<IMS,SPB>
{
  using Base = DynamicPowerPreBasis<IMS,SPB>;

public:

  //! The child pre-basis
  using SubPreBasis = SPB;

  //! Template mapping root tree path to type of created tree node
  using Node = PowerBasisNode<typename SubPreBasis::Node, C>;

  //! Type used for indices and size information
  using size_type = typename Base::size_type;

  //! Strategy used to merge the global indices of the child factories
  using IndexMergingStrategy = IMS;

  //! Number of children provided as an integral constant
  inline static constexpr std::integral_constant<std::size_t,C> children = {};

  /**
   * \brief Constructor for given child pre-basis objects for static size of the power-basis
   *
   * The child factories will be stored as copies
   */
  template<class... SFArgs,
    disableCopyMove<PowerPreBasis, SFArgs...> = 0,
    enableIfConstructible<SubPreBasis, SFArgs...> = 0>
  explicit PowerPreBasis(SFArgs&&... sfArgs) :
    Base(std::size_t(C), std::forward<SFArgs>(sfArgs)...)
  {}

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    Node node{};
    for (std::size_t i=0; i<children(); ++i)
      node.setChild(i, Base::subPreBasis_.makeNode());
    return node;
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return size(Dune::ReservedVector<size_type, Base::multiIndexBufferSize>{});
  }

  //! Return number of possible values for next position in multi index
  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
  {
    return Base::sizeImpl(prefix, children, IndexMergingStrategy{});
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<class NodeType, typename It,
    std::enable_if_t<NodeType::isPower, int> = 0>
  It indices(const NodeType& node, It it) const
  {
    return Base::indicesImpl(node, it, children, IndexMergingStrategy{});
  }

  //! Return the associated container descriptor
  auto containerDescriptor() const
  {
    return Base::containerDescriptorImpl(children);
  }
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
[[deprecated("Using the method `power` without an explicit index merging strategy"
  " will change its meaning after the release of dune-functions 2.11."
  " Previously, the default merging strategy was `BlockedInterleaved`,"
  " but this will change to `FlatInterleaved`."
  " Change the call to `power(..., blockedInterleaved())` to retain the old behavior.")]]
auto power(ChildPreBasisFactory&& childPreBasisFactory)
{
  return [childPreBasisFactory](const auto& gridView) {
    auto childPreBasis = childPreBasisFactory(gridView);
    return PowerPreBasis<BlockedInterleaved, decltype(childPreBasis), k>(std::move(childPreBasis));
  };
}

} // end namespace BasisFactory


} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_POWERBASIS_HH
