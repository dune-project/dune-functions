// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_GEOMETRYINANCESTOR_HH
#define DUNE_FUNCTIONS_COMMON_GEOMETRYINANCESTOR_HH

#include <vector>

#include <dune/common/iteratorrange.hh>

namespace Dune::Functions {



/**
 * \brief A geometry embedding a descendent element into an ancestor
 *
 * \ingroup FunctionUtility
 *
 * \tparam Element Type of elements considered by this embedding
 *
 * This class will store a chain of `geometryInFather` objects
 * connecting a coarse element with a fine element. The `global()`
 * and `local()` methods of the `GeometryInAncestor` simply chain
 * those stored `geometryInFather` objects.
 * Since this requires that objects of type `GeometryInAncestor` allocate
 * dynamic memory, they are not intended to be created on the fly
 * within a hot loop. Instead one reused a single `GeometryInAncestor`
 * object by binding if to a new element.
 *
 * Currently this only provides the `local()` and `global()` methods
 * of the geometry interface.
 */
template<class Element>
class GeometryInAncestor
{
  using GeometryInFather = typename Element::LocalGeometry;

  template<class V>
  static auto reverse_view(V&& v)
  {
    return Dune::IteratorRange(v.rbegin(), v.rend());
  }

public:

  //! Type of local coordinate (local within fine element)
  using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
  //
  //! Type of global coordinate (local within coarse element)
  using GlobalCoordinate = typename Element::Geometry::LocalCoordinate;

  GeometryInAncestor() = default;
  GeometryInAncestor(GeometryInAncestor&&) = delete;

  /**
   * \brief Copy constructor
   *
   * If the original GeometryInAncestor holds a pointer
   * to an external fine Element, the new GeometryInAncestor
   * points to the same external fine Element.
   */
  GeometryInAncestor(const GeometryInAncestor& other)
    : storedFineElement_(other.storedFineElement_)
    , storedCoarseElement_(other.storedCoarseElement_)
    , fineElementPtr_(&storedFineElement_)
    , coarseElementPtr_(&storedFineElement_)
    , geometryInFathersVector_(other.geometryInFathersVector_)
  {
    // If other points to external fine element,
    // we point to the same one.
    if (other.fineElementPtr_ != & other.storedFineElement_)
      fineElementPtr_ = other.fineElementPtr_;
    // If the other coarse element is not the stored one,
    // then it's the fine one.
    if (other.coarseElementPtr_ != & other.storedCoarseElement_)
      coarseElementPtr_ = fineElementPtr_;
  }

  /**
   * \brief Copy construction setting an external fine element
   *
   * This copies the original GeometryInAncestor and sets the
   * pointer to the external fine element passed additionally.
   */
  GeometryInAncestor(const GeometryInAncestor& other, const Element& fineElement)
    : storedFineElement_(other.storedFineElement_)
    , storedCoarseElement_(other.storedCoarseElement_)
    , fineElementPtr_(&fineElement)
    , coarseElementPtr_(&storedFineElement_)
    , geometryInFathersVector_(other.geometryInFathersVector_)
  {
    // If the other coarse element is not the stored one,
    // then it's the fine one.
    if (other.coarseElementPtr_ != & other.storedCoarseElement_)
      coarseElementPtr_ = fineElementPtr_;
  }

  /**
   * \brief Bind this GeometryInAncestor to a fine element and search coarse element
   *
   * \param fineElement The fine element to be bound to
   * \param includeFather Unary predicate indicating whether father element should be traversed
   *
   * Build the geometry information by traversing the fathers starting from
   * the given `fineElement` as long as `includeFather(element)` returns `true`
   * to obtain a chain of `geometryInFather` objects connecting the `fineElement`
   * with a `coarseElement`. The `coarseElement` is the first element where
   * `includeFather(coarseElement)` returns `false`.
   *
   * This overload for an l-value `fineElement` will store
   * a pointer to `fineElement`.
   */
  template<class F>
  requires (std::is_invocable_r_v<bool,F,Element>)
  const Element& bind(const Element& fineElement, F&& includeFather)
  {
    fineElementPtr_ = &fineElement;
    geometryInFathersVector_.clear();
    coarseElementPtr_ = fineElementPtr_;
    while (includeFather(*coarseElementPtr_))
    {
      geometryInFathersVector_.push_back(coarseElementPtr_->geometryInFather());
      storedCoarseElement_ = coarseElementPtr_->father();
      coarseElementPtr_ = &storedCoarseElement_;
    }
    return coarseElement();
  }

  /**
   * \brief Bind this GeometryInAncestor to a fine element and search coarse element
   *
   * \param fineElement The fine element to be bound to
   * \param includeFather Unary predicate indicating whether father element should be traversed
   *
   * Build the geometry information by traversing the fathers starting from
   * the given `fineElement` as long as `includeFather(element)` returns `true`
   * to obtain a chain of `geometryInFather` objects connecting the `fineElement`
   * with a `coarseElement`. The `coarseElement` is the first element where
   * `includeFather(coarseElement)` returns `false`.
   *
   * This overload for an r-value `fineElement` will store
   * a copy of the `fineElement`.
   */
  template<class F>
  requires (std::is_invocable_r_v<bool,F,Element>)
  const Element& bind(Element&& fineElement, F&& includeFather)
  {
    storedFineElement_ = fineElement;
    bind(storedFineElement_, includeFather);
  }

  //! Return the fine element
  const Element& fineElement() const
  {
    return *fineElementPtr_;
  }

  //! Return the coarse element
  const Element& coarseElement() const
  {
    return *coarseElementPtr_;
  }

  //! Map local coordinate in fine element into coarse element
  GlobalCoordinate global(LocalCoordinate x) const
  {
    for (const auto& g : geometryInFathersVector_)
      x = g.global(x);
    return x;
  }

  //! Map local coordinate in coarse element into fine element
  LocalCoordinate local(GlobalCoordinate x) const
  {
    // Here we can be instead use std::ranges::reverse_view once
    // all supported compilers provide this (gcc>=10, clang>=16).
    for (const auto& g : reverse_view(geometryInFathersVector_))
      x = g.local(x);
    return x;
  }

private:
  Element storedFineElement_;
  Element storedCoarseElement_;
  const Element* fineElementPtr_ = nullptr;
  const Element* coarseElementPtr_ = nullptr;
  std::vector<GeometryInFather> geometryInFathersVector_;
};



} // namespace Dune::Functions




#endif // DUNE_FUNCTIONS_COMMON_GEOMETRYINANCESTOR_HH
