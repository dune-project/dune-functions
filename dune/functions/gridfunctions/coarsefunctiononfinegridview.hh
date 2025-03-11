// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_COARSEFUNCTIONONFINEGRIDVIEW_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_COARSEFUNCTIONONFINEGRIDVIEW_HH

#include <optional>
#include <type_traits>
#include <utility>

#include <dune/common/referencehelper.hh>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/geometryinancestor.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune::Functions {



/**
 * \brief A wrapper representing a coarse grid function on a fine gridview
 *
 * \ingroup FunctionImplementations
 *
 * \tparam GridFunction Type of the wrapped grid function
 * \tparam GV Type of the target grid view this function should act on
 *
 * This wraps a grid function such that it can be used as a `GridViewFunction`
 * on a user-provided `GridView` under the following assumptions:
 * 1. The grid function's entity set and the `GridView` belong to the same grid.
 * 2. The entity set is coarser than the `GridView` in the sense that any
 *    element from the `GridView` has an ancestor in the entity set.
 */
template<class GridFunction, class GV, template<class> class DerivativeTraits=Dune::Functions::DefaultDerivativeTraits>
class CoarseFunctionOnFineGridView
{
  using RawGridFunction = Dune::ResolveRef_t<GridFunction>;

  auto&& rawFunction() const
  {
    return Dune::resolveRef(function_);
  }

public:

  using GridView = GV;
  using EntitySet = Dune::Functions::GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;
  using Domain = typename EntitySet::GlobalCoordinate;
  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Range = std::decay_t<decltype(std::declval<RawGridFunction>()(std::declval<Domain>()))>;

private:

  using CoarseEntitySet = std::decay_t<decltype(std::declval<RawGridFunction>().entitySet())>;
  using GeometryInAncestor = Dune::Functions::GeometryInAncestor<Element>;
  using Traits = Dune::Functions::Imp::GridFunctionTraits<Range(Domain), EntitySet, DerivativeTraits, 56>;

  class CoarseLocalFunctionOnFineGridView
  {
    using Traits = typename CoarseFunctionOnFineGridView::Traits::LocalFunctionTraits;

  public:

    using Derivative = decltype(localFunction(derivative(std::declval<CoarseFunctionOnFineGridView>())));
    using RawLocalFunction = std::decay_t<decltype(localFunction(std::declval<const RawGridFunction&>()))>;

    /**
     * \brief Construct the LocalFunction
     *
     * The LocalFunction is created from the global CoarseFunctionOnFineGridView.
     */
    CoarseLocalFunctionOnFineGridView(RawLocalFunction&& localFunction, const CoarseEntitySet& coarseEntitySet)
      : element_()
      , localFunction_(localFunction)
      , coarseEntitySet_(coarseEntitySet)
      , geometryInAncestor_()
    {}

    /**
     * \brief Construct the LocalFunction
     *
     * The LocalFunction is created from the global CoarseFunctionOnFineGridView.
     */
    CoarseLocalFunctionOnFineGridView(
        RawLocalFunction&& localFunction,
        const CoarseEntitySet& coarseEntitySet,
        const GeometryInAncestor& geometryInAncestor,
        const std::optional<Element>& element
      )
      : element_(element)
      , localFunction_(localFunction)
      , coarseEntitySet_(coarseEntitySet)
      , geometryInAncestor_(geometryInAncestor, *element_)
    {}

    //! Bind to an element from the GridView
    void bind(const Element& element)
    {
      element_ = element;
      geometryInAncestor_.bind(*element_, [&](const auto& e) { return not coarseEntitySet_.contains(e); });
      localFunction_.bind(geometryInAncestor_.coarseElement());
    }

    //! \brief Unbind
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
    friend auto derivative(const CoarseLocalFunctionOnFineGridView& f)
    {
      if constexpr(requires{ derivative(f.localFunction_); })
        return Derivative(derivative(f.localFunction_), f.coarseEntitySet_, f.geometryInAncestor_, f.element_);
      else
        return typename Traits::DerivativeInterface{};
    }

    //! Evaluate function in local coordinates
    Range operator()(LocalDomain x) const
    {
      return localFunction_(geometryInAncestor_.global(x));
    }

  private:
    std::optional<Element> element_;
    RawLocalFunction localFunction_;
    const CoarseEntitySet& coarseEntitySet_;
    GeometryInAncestor geometryInAncestor_;
  };

public:

  using LocalFunction = CoarseLocalFunctionOnFineGridView;

  /**
   * \brief Create CoarseFunctionOnFineGridView from GridFunction and GridView
   *
   * \param gridFunction The GridFunction that should be represented on gridView
   * \param gridView The GridFunction should be represented on this gridView
   */
  CoarseFunctionOnFineGridView(const GridFunction& function, const GridView& gridView)
    : function_(function)
    , entitySet_(gridView)
  {}

  /**
   * \brief Create CoarseFunctionOnFineGridView from GridFunction and GridView
   *
   * \param gridFunction The GridFunction that should be represented on gridView
   * \param gridView The GridFunction should be represented on this gridView
   */
  CoarseFunctionOnFineGridView(GridFunction&& function, const GridView& gridView)
    : function_(std::move(function))
    , entitySet_(gridView)
  {}

  //! Evaluate function in global coordinates
  Range operator()(const Domain& x) const
  {
    return function_(x);
  }

  //! Obtain global derivative of this function
  friend auto derivative(const CoarseFunctionOnFineGridView& f)
  {
    if constexpr(requires{ derivative(f.rawFunction()); })
    {
      using RawDerivative = std::decay_t<decltype(derivative(f.rawFunction()))>;
      return CoarseFunctionOnFineGridView<RawDerivative, GridView, DerivativeTraits>(derivative(f.rawFunction()), f.entitySet_.gridView());
    }
    else
      return typename Traits::DerivativeInterface{};
  }

  //! Create a LocalFunction for evaluation in local coordinates
  friend LocalFunction localFunction(const CoarseFunctionOnFineGridView& f)
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

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_COARSEFUNCTIONONFINEGRIDVIEW_HH
