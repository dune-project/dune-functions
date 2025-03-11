// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_FINEFUNCTIONONCOARSEGRIDVIEW_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_FINEFUNCTIONONCOARSEGRIDVIEW_HH

#include <optional>
#include <type_traits>
#include <utility>
#include <limits>
#include <algorithm>
#include <cmath>

#include <dune/common/referencehelper.hh>

#include <dune/geometry/type.hh>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune::Functions {

namespace Impl {

namespace ReferenceElementUtilities {

// Compute the l1-distance of x to the reference element identified
// by the topology id and dimension.
template<class X, class FT = std::decay_t<decltype(std::declval<X>()[0])> >
FT distance(unsigned int topologyId, int dim, X x, FT scaleFactor = FT(1))
{
  using std::abs;
  using std::max;
  using std::min;
  auto dist_x_last = max(max(x[dim-1]-scaleFactor, -x[dim-1]), FT(0));
  if (dim > 1)
  {
    if (Dune::Impl::isPyramid(topologyId, dim))
      scaleFactor -= max(min(x[dim-1], scaleFactor), FT(0));
    return distance(Dune::Impl::baseTopologyId(topologyId, dim), dim-1, x, scaleFactor) + dist_x_last;
  }
  if (dim == 1)
    return dist_x_last;
  return FT(0);
}

// Check if the l1-distance of x to the reference element identified
// by the topology id and dimension is less than a tolerance. This
// implementation is significantly faster than checking if
// distance(...) <= tolerance. It is almost as fast as checkInside(...)
// of the refenece element, but the latter does not reflect the
// distance wrt any norm while we use the l1-norm here.
template<class X, class FT = std::decay_t<decltype(std::declval<X>()[0])> >
bool checkInside(unsigned int topologyId, int dim, X x, FT tolerance, FT scaleFactor = FT(1))
{
  using std::abs;
  using std::max;
  using std::min;
  if (dim > 0)
  {
    auto dist_x_last = max(x[dim-1]-scaleFactor, -x[dim-1]);
    if (dist_x_last <= tolerance)
    {
      if (Dune::Impl::isPyramid(topologyId, dim))
        scaleFactor -= max(min(x[dim-1], scaleFactor), FT(0));
      return checkInside(Dune::Impl::baseTopologyId(topologyId, dim), dim-1, x, tolerance - max(dist_x_last, FT(0)), scaleFactor);
    }
    return false;
  }
  return true;
}

} // namespace ReferenceElementUtilities

} // namespace Impl







/**
 * \brief A wrapper representing a fine grid function on a gridview
 *
 * \ingroup FunctionImplementations
 *
 * \tparam GridFunction Type of the wrapped grid function
 * \tparam GV Type of the target grid view this function should act on
 *
 * This wraps a grid function such that it can be used as a `GridViewFunction`
 * on a user-provided `GridView` under the following assumptions:
 * 1. The grid function's entity set and the `GridView` belong to the same grid.
 * 2. The entity set is finer than the `GridView` in the sense that any
 *    element from the entity has an ancestor in the `GridView`.
 */
template<class GridFunction, class GV, template<class> class DerivativeTraits=Dune::Functions::DefaultDerivativeTraits>
class FineFunctionOnCoarseGridView
{
  using RawGridFunction = Dune::ResolveRef_t<GridFunction>;

  auto&& rawFunction() const
  {
    return Dune::resolveRef(function_);
  }

  static constexpr auto dim = GV::Grid::dimension;

public:

  using GridView = GV;
  using EntitySet = Dune::Functions::GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;
  using Domain = typename EntitySet::GlobalCoordinate;
  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Range = std::decay_t<decltype(std::declval<RawGridFunction>()(std::declval<Domain>()))>;

private:

  using FineEntitySet = std::decay_t<decltype(std::declval<RawGridFunction>().entitySet())>;
  using Traits = Dune::Functions::Imp::GridFunctionTraits<Range(Domain), EntitySet, DerivativeTraits, 56>;

  class FineLocalFunctionOnCoarseGridView
  {
    using Traits = typename FineFunctionOnCoarseGridView::Traits::LocalFunctionTraits;

  public:

    using Derivative = decltype(localFunction(derivative(std::declval<FineFunctionOnCoarseGridView>())));
    using RawLocalFunction = std::decay_t<decltype(localFunction(std::declval<const RawGridFunction&>()))>;

    /**
     * \brief Construct the LocalFunction
     *
     * The LocalFunction is created from the global FineFunctionOnCoarseGridView.
     **/
    FineLocalFunctionOnCoarseGridView(RawLocalFunction&& localFunction, const FineEntitySet& fineEntitySet)
      : element_()
      , localFunction_(localFunction)
      , fineEntitySet_(fineEntitySet)
      , forwardToFineFunction_(false)
    {}

    /**
     * \brief Construct the LocalFunction
     *
     * The LocalFunction is created from the global FineFunctionOnCoarseGridView.
     **/
    FineLocalFunctionOnCoarseGridView(
        RawLocalFunction&& localFunction,
        const FineEntitySet& fineEntitySet,
        bool forwardToFineFunction,
        const std::optional<Element>& element
      )
      : element_(element)
      , localFunction_(localFunction)
      , fineEntitySet_(fineEntitySet)
      , forwardToFineFunction_(forwardToFineFunction)
    {}

    //! Bind to an element from the GridView
    void bind(const Element& element)
    {
      element_ = element;
      forwardToFineFunction_ = fineEntitySet_.contains(*element_);
      if (forwardToFineFunction_)
        localFunction_.bind(element);
    }

    //! \brief Unbind the inner local-functions.
    void unbind()
    {
      element_.reset();
    }

    //! Return if the local function is bound to an element of the GridView
    bool bound() const
    {
      return static_cast<bool>(element_);
    }

    //! Obtain the grid element this function is bound to
    const Element& localContext() const
    {
      return *element_;
    }

    //! Obtain local derivative of this function
    friend auto derivative(const FineLocalFunctionOnCoarseGridView& f)
    {
      if constexpr(requires{ derivative(f.localFunction_); })
        return Derivative(derivative(f.localFunction_), f.fineEntitySet_, f.forwardToFineFunction_, f.element_);
      else
        return typename Traits::DerivativeInterface{};
    }

    //! Evaluate function in local coordinates
    Range operator()(LocalDomain x) const
    {
      if (forwardToFineFunction_)
        return localFunction_(x);
      return evaluateInDescendent(*element_, x);
    }

  private:

    // Find a child containing the point and evaluate there recursively
    Range evaluateInDescendent(const Element& element, LocalDomain x) const
    {
      Element closestChild;
      LocalDomain xInClosestChild;
      double distanceToClosestChild = std::numeric_limits<double>::max();
      for(const auto& child : descendantElements(element, element.level()+1))
      {
        auto&& geometry = child.geometryInFather();
        auto xInChild = geometry.local(x);
        auto dist = Impl::ReferenceElementUtilities::distance(child.type().id(), dim, xInChild);
        if (dist < distanceToClosestChild)
        {
          closestChild = child;
          distanceToClosestChild = dist;
          xInClosestChild = xInChild;
          if (distanceToClosestChild==0)
            break;
        }
      }
      if (fineEntitySet_.contains(closestChild))
      {
        localFunction_.bind(closestChild);
        return localFunction_(xInClosestChild);
      }
      else
        return evaluateInDescendent(closestChild, xInClosestChild);
    }

    std::optional<Element> element_;
    mutable RawLocalFunction localFunction_;
    const FineEntitySet& fineEntitySet_;
    bool forwardToFineFunction_ = false;
  };

public:

  using LocalFunction = FineLocalFunctionOnCoarseGridView;

  /**
   * \brief Create FineFunctionOnCoarseGridView from GridFunction and GridView
   *
   * \param gridFunction The GridFunction that should be represented on gridView
   * \param gridView The GridFunction should be represented on this gridView
   */
  FineFunctionOnCoarseGridView(const GridFunction& function, const GridView& gridView)
    : function_(function)
    , entitySet_(gridView)
  {}

  /**
   * \brief Create FineFunctionOnCoarseGridView from GridFunction and GridView
   *
   * \param gridFunction The GridFunction that should be represented on gridView
   * \param gridView The GridFunction should be represented on this gridView
   */
  FineFunctionOnCoarseGridView(GridFunction&& function, const GridView& gridView)
    : function_(std::move(function))
    , entitySet_(gridView)
  {}

  //! Evaluate function in global coordinates
  Range operator()(const Domain& x) const
  {
    return function_(x);
  }

  //! Obtain global derivative of this function
  friend auto derivative(const FineFunctionOnCoarseGridView& f)
  {
    if constexpr(requires{ derivative(f.rawFunction()); })
    {
      using RawDerivative = std::decay_t<decltype(derivative(f.rawFunction()))>;
      return FineFunctionOnCoarseGridView<RawDerivative, GridView, DerivativeTraits>(derivative(f.rawFunction()), f.entitySet_.gridView());
    }
    else
      return typename Traits::DerivativeInterface{};
  }

  //! Create a LocalFunction for evaluation in local coordinates
  friend LocalFunction localFunction(const FineFunctionOnCoarseGridView& f)
  {
    return LocalFunction(localFunction(f.rawFunction()), f.rawFunction().entitySet());
  }

  //! Return the EntitySet associated to this GridViewFunction
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

protected:

  GridFunction function_;
  EntitySet entitySet_;
};



} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_FINEFUNCTIONONCOARSEGRIDVIEW_HH
