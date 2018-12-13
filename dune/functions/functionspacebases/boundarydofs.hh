// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BOUNDARYDOFS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BOUNDARYDOFS_HH

#include <utility>

#include <dune/functions/functionspacebases/subentitydofs.hh>

namespace Dune {
namespace Functions {



/**
 * \brief Loop over all DOFs on the boundary
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This loops over all DOFs of a basis associated to sub-entities
 * on the boundary. This overload will pass three arguments to the
 * given loop callback: The local index of the boundary DOF,
 * a bound local view this local index belongs to,
 * and a boundary intersection associated to a sub-entity such that
 * the DOF is associated to a sub-sub-entity of this sub-entity.
 * Notice that this may visit the same DOF multiple times.
 *
 * If this callback signature is not suitable you can use one of
 * the another variants of forEachBoundaryDOF.
 *
 * \param basis A function space basis
 * \param f A callback that will be called with a local index, a bound local view, and an intersection of the visited boundary DOF
 */
template<class Basis, class F,
  decltype(std::declval<std::decay_t<F>>()(0, std::declval<typename Basis::LocalView>(),std::declval<typename Basis::GridView::Intersection>()), 0) = 0>
void forEachBoundaryDOF(const Basis& basis, F&& f)
{
  auto localView = basis.localView();
  auto seDOFs = subEntityDOFs(basis);
  const auto& gridView = basis.gridView();
  for(auto&& element : elements(gridView))
    if (element.hasBoundaryIntersections())
    {
      localView.bind(element);
      for(const auto& intersection: intersections(gridView, element))
        if (intersection.boundary())
          for(auto localIndex: seDOFs.bind(localView,intersection))
            f(localIndex, localView, intersection);
    }
}



/**
 * \brief Loop over all DOFs on the boundary
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This loops over all DOFs of a basis associated to sub-entities
 * on the boundary. This overload will pass two arguments to the
 * given loop callback: The local index of the boundary DOF and
 * a bound local view this local index belongs to.
 * Notice that this may visit the same DOF multiple times.
 *
 * If this callback signature is not suitable you can use one of
 * the another variants of forEachBoundaryDOF.
 *
 * \param basis A function space basis
 * \param f A callback that will be called with a local index and a bound local view of the visited boundary DOF
 */
template<class Basis, class F,
  decltype(std::declval<std::decay_t<F>>()(0, std::declval<typename Basis::LocalView>()),0) = 0>
void forEachBoundaryDOF(const Basis& basis, F&& f)
{
  auto localView = basis.localView();
  auto seDOFs = subEntityDOFs(basis);
  const auto& gridView = basis.gridView();
  for(auto&& element : elements(gridView))
    if (element.hasBoundaryIntersections())
    {
      localView.bind(element);
      for(const auto& intersection: intersections(gridView, element))
        if (intersection.boundary())
          for(auto localIndex: seDOFs.bind(localView,intersection))
            f(localIndex, localView);
    }
}



/**
 * \brief Loop over all DOFs on the boundary
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This loops over all DOFs of a basis associated to sub-entities
 * on the boundary. This overload will pass a single arguments to the
 * given loop callback: The global (multi-)index of the boundary DOF.
 * Notice that this may visit the same DOF multiple times.
 *
 * If this callback signature is not suitable you can use one of
 * the another variants of forEachBoundaryDOF.
 *
 * \param basis A function space basis
 * \param f A callback that will be called with the global index of the visited boundary DOF
 */
template<class Basis, class F,
  decltype(std::declval<std::decay_t<F>>()(std::declval<typename Basis::MultiIndex>()),0) = 0>
void forEachBoundaryDOF(const Basis& basis, F&& f)
{
  auto localView = basis.localView();
  auto seDOFs = subEntityDOFs(basis);
  const auto& gridView = basis.gridView();
  for(auto&& element : elements(gridView))
    if (element.hasBoundaryIntersections())
    {
      localView.bind(element);
      for(const auto& intersection: intersections(gridView, element))
        if (intersection.boundary())
          for(auto localIndex: seDOFs.bind(localView,intersection))
            f(localView.index(localIndex));
    }
}



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BOUNDARYDOFS_HH
