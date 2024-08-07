// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH

#include <memory>

#include <dune/common/concept.hh>

#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>


namespace Dune {
namespace Functions {



template<class Signature, class GridView, template<class> class DerivativeTraits=DefaultDerivativeTraits, size_t bufferSize=56>
class GridViewFunction
{};



/**
 * \brief Wrapper class for functions defined on a GridView.
 *
 * \ingroup FunctionInterface
 *
 * Being defined on a grid view means in particular that you can evaluate the function
 * in local coordinates of a given element of the grid view.
 *
 * This models the \ref Concept::GridViewFunction<Range(Domain), GridView, DerivativeTraits> concept.
 *
 * \tparam GV The GridView that the function is defined on
 * \tparam Domain The domain type used for function arguments
 * \tparam Range The range type used for function values
 */
template<class Range, class Domain, class GV, template<class> class DerivativeTraits, size_t bufferSize>
class GridViewFunction<Range(Domain), GV, DerivativeTraits, bufferSize> :
  public GridFunction<Range(Domain), GridViewEntitySet<GV, 0>, DerivativeTraits, bufferSize>
{
  using Base = GridFunction<Range(Domain), GridViewEntitySet<GV, 0>, DerivativeTraits, bufferSize>;
public:
  using GridView = GV;

  using Base::Base;
};



/**
 * \brief Construct a function modeling GridViewFunction from function and grid view.
 *
 * This specialization is used for functions that already
 * support `localFunction()`. It will simply return a copy of `f`.
 *
 * \param f A function object supporting argument compatible with global coordinates
 * \param gridView The GridView the function should act on.
 *
 * \returns A function that models the GridViewFunction interface.
 */
template<class F, class GridView,
  std::enable_if_t<
    models< Imp::HasFreeLocalFunction, F>() , int> = 0>
std::decay_t<F>
  makeGridViewFunction(F&& f, const GridView& gridView)
{
  return std::forward<F>(f);
}



/**
 * \brief Construct a function modeling GridViewFunction from function and grid view.
 *
 * This specialization is used for functions that do not
 * support `localFunction()` themselves. It will forward
 * to `makeAnalyticGridViewFunction()`.
 *
 * Notice that the returned function will store a copy of
 * the original function and the GridView.
 *
 * \param f A function object supporting argument compatible with global coordinates
 * \param gridView The GridView the function should act on.
 *
 * \returns A function that models the GridFunction interface.
 */
template<class F, class GridView,
  std::enable_if_t<
    not(models< Imp::HasFreeLocalFunction, F>()) , int> = 0>
auto makeGridViewFunction(F&& f, const GridView& gridView)
  -> decltype(makeAnalyticGridViewFunction(std::forward<F>(f), gridView))
{
  return makeAnalyticGridViewFunction(std::forward<F>(f), gridView);
}



} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH
