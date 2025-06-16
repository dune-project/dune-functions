// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH


#include <tuple>
#include <optional>

#include <dune/common/concept.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>

#include <dune/functions/common/overflowarray.hh>
#include <dune/functions/common/multiindex.hh>
#include <dune/functions/functionspacebases/concepts.hh>



namespace Dune {
namespace Functions {



/** \brief The restriction of a finite element basis to a single element */
template<class GB>
class DefaultLocalView
{
public:

  //! The global FE basis that this is a view on
  using GlobalBasis = GB;

  //! The grid view the global FE basis lives on
  using GridView = typename GlobalBasis::GridView;

  //! Type of the grid element we are bound to
  using Element = typename GridView::template Codim<0>::Entity;

  //! The type used for sizes
  using size_type = std::size_t;

  //! Tree of local finite elements / local shape function sets
  using Tree = typename GlobalBasis::PreBasis::Node;

protected:

  using PreBasis = typename GlobalBasis::PreBasis;

  // Type used to store the multi indices of the basis vectors.
  // In contrast to MultiIndex this always has dynamic size.
  // It's guaranteed, that you can always cast it to MultiIndex
  using MultiIndexStorage =
      std::conditional_t<(PreBasis::minMultiIndexSize == PreBasis::maxMultiIndexSize),
        OverflowArray<StaticMultiIndex<size_type, PreBasis::maxMultiIndexSize>, PreBasis::multiIndexBufferSize>,
        Dune::ReservedVector<size_type, PreBasis::multiIndexBufferSize>>;

public:

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex =
      std::conditional_t<(PreBasis::minMultiIndexSize == PreBasis::maxMultiIndexSize),
        StaticMultiIndex<size_type, PreBasis::maxMultiIndexSize>,
        Dune::ReservedVector<size_type, PreBasis::multiIndexBufferSize>>;


  /** \brief Construct local view for a given global finite element basis */
  DefaultLocalView(const GlobalBasis& globalBasis) :
    globalBasis_(&globalBasis),
    tree_(globalBasis_->preBasis().makeNode())
  {
    static_assert(models<Concept::BasisTree<GridView>, Tree>(), "Tree type passed to DefaultLocalView does not model the BasisNode concept.");
    initializeTree(tree_);
  }

  /** \brief Deep copy of the local view */
  DefaultLocalView(const DefaultLocalView& other) :
    globalBasis_(other.globalBasis_),
    tree_(globalBasis_->preBasis().makeNode())
  {
    initializeTree(tree_);
    if (other.bound())
      bind(other.element());
  }

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = e;
    bindTree(tree_, *element_);
    indices_.resize(size());
    globalBasis_->preBasis().indices(tree_, indices_.begin());
  }

  /** \brief Return if the view is bound to a grid element
   */
  bool bound() const
  {
    return static_cast<bool>(element_);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   */
  void unbind()
  {
    element_.reset();
  }

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /** \brief Total number of degrees of freedom on this element
   */
  size_type size() const
  {
    return tree_.size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   */
  size_type maxSize() const
  {
    return globalBasis_->preBasis().maxNodeSize();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  const MultiIndex& index(size_type i) const
  {
    return indices_[i];
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

  const DefaultLocalView& rootLocalView() const
  {
    return *this;
  }

protected:
  const GlobalBasis* globalBasis_;
  std::optional<Element> element_;
  Tree tree_;
  std::vector<MultiIndexStorage> indices_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH
