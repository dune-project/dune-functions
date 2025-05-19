// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RESTRICTEDBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RESTRICTEDBASIS_HH

#include <cstddef>
#include <utility>

#include <dune/typetree/traversal.hh>

#include <dune/functions/functionspacebases/nodes.hh>


namespace Dune {
namespace Functions {



/**
 * \brief A pre-basis restricted to a sub-domain
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * This pre-basis wraps another pre-basis and restricts it to a sub-domain.
 * The wrapped pre-basis is assumed to be defined on a sub-domain of the
 * full grid view. Then the RestrictedPreBasis defines a pre-basis on
 * the full grid view implemented in terms of the pre-basis on the sub-domain.
 * Most of the methods are forwarded to the sub-domain pre-basis with
 * two exceptions: When binding a `RestrictedPreBasis::Node`, this only
 * calls `subDomainNode.bind(element)` if the element is contained in
 * the sub-domain, otherwise the sizes of the node and all its descendents
 * is set to zero. Furthermore it only calls `subDomainPreBasis.indices(node,it)`
 * if `node.size()` is not zero.
 *
 * \tparam GV Type of the (full) grid view this pre-basis is defined on.
 * \tparam SDPB Type of a pre-basis defined on sub-set of the full grid view
 * \tparam SD Type of the sub-domain
 */
template<class GV, class SDPB, class SD>
class RestrictedPreBasis
{
  using This = RestrictedPreBasis<GV, SDPB, SD>;

public:

  using SubDomain = SD;
  using SubDomainPreBasis = SDPB;
  using SubDomainGridView = typename SubDomainPreBasis::GridView;

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  class Node
    : public SubDomainPreBasis::Node
  {
    using Base = typename SubDomainPreBasis::Node;
  public:

    using Element = typename Base::Element;

    Node(Base&& base, const SubDomainGridView& subDomainGridView, const SubDomain& subDomain)
      : Base(base)
      , subDomainGridView_(subDomainGridView)
      , subDomain_(subDomain)
    {}

    void bind(const Element& element)
    {
      if (subDomainGridView_.contains(element))
        Base::bind(element);
      else
      {
        Dune::TypeTree::forEachNode(static_cast<Base&>(*this) , [&](auto& node, const auto& treePath) {
          Dune::Functions::Impl::BasisNodeSetupHelper::setOffset(node, this->offset());
          Dune::Functions::Impl::BasisNodeSetupHelper::setSize(node, 0);
        });
      }
    }

    const SubDomain& subDomain() const
    {
      return subDomain_;
    }

  private:
    const SubDomainGridView& subDomainGridView_;
    const SubDomain& subDomain_;
  };

  static constexpr size_type maxMultiIndexSize = SubDomainPreBasis::maxMultiIndexSize;
  static constexpr size_type minMultiIndexSize = SubDomainPreBasis::minMultiIndexSize;
  static constexpr size_type multiIndexBufferSize = SubDomainPreBasis::multiIndexBufferSize;

  /**
   * \brief Constructor for given sub-domain pre-basis
   *
   * The grid view and sub-domain pre-basis will be stored as copy
   * while a pointer to the sub-domain object is stored.
   */
  RestrictedPreBasis(const GridView& gridView, SubDomainPreBasis&& subDomainPreBasis, const SubDomain& subDomain)
    : gridView_(gridView)
    , subDomainPreBasis_(std::move(subDomainPreBasis))
    , subDomainPtr_(&subDomain)
  {}

  //! Initialize the global indices
  void initializeIndices()
  {
    subDomainPreBasis_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return gridView_;
  }

  /**
   * \brief Update the stored grid view, to be called
   *
   * This will also call `subDomainPreBasis.update(subDomainPtr->gridView())`
   * with the stored sub-domain pointer to update the sub-domain pre-basis.
   * Hence it requires the `subDomain` object passed to the present pre-basis
   * has been correctly updated externally, before calling this method.
   */
  void update(const GridView& gv)
  {
    gridView_ = gv;
    subDomainPreBasis_.update(subDomainPtr_->gridView());
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
    return Node(subDomainPreBasis_.makeNode(), subDomainPreBasis_.gridView(), *subDomainPtr_);
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return subDomainPreBasis_.size();
  }

  //! Return number of possible values for next position in multi index
  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
  {
    return subDomainPreBasis_.size(prefix);
  }

  //! Return the container descriptor of the pre-basis
  auto containerDescriptor() const
  {
    return subDomainPreBasis_.containerDescriptor();
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return subDomainPreBasis_.dimension();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return subDomainPreBasis_.maxNodeSize();
  }

  const SubDomainPreBasis& subDomainPreBasis() const
  {
    return subDomainPreBasis_;
  }

  SubDomainPreBasis& subDomainPreBasis()
  {
    return subDomainPreBasis_;
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    if (node.size() == 0)
      return it;
    else
      return subDomainPreBasis_.indices(node, it);
  }

protected:
  GridView gridView_;
  SubDomainPreBasis subDomainPreBasis_;
  const SubDomain* subDomainPtr_;
};



namespace BasisFactory {

/**
 * \brief Create a RestrictedPreBasisFactory
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * This creates a preBasisFactory that restricts a pre-basis
 * to a sub-domain. The \p subDomain object must support
 * `subDomain.gridView()`. The returned `gridView` must
 * be defined on a subset of the host grid view, it must
 * support `gridView.contains(element)` to check if an
 * element is contained in the sub-domain and it must provide
 * all the interface functions required by the to-be-constructed
 * pre-basis type.
 *
 * \param subPreBasisFactory A PreBasisFactory use to create a pre-basis on the sub-domain
 * \param subDomain A sub-domain object
 */
template<class SubDomainBasisFactory, class SubDomain>
auto restrict(SubDomainBasisFactory&& subPreBasisFactory, const SubDomain& subDomain)
{
  return [
    subPreBasisFactory=std::forward<SubDomainBasisFactory>(subPreBasisFactory),
    &subDomain
  ](const auto& gridView) {
    return Dune::Functions::RestrictedPreBasis(gridView, subPreBasisFactory(subDomain.gridView()), subDomain);
  };
}


} // end namespace BasisFactory
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RESTRICTEDBASIS_HH
