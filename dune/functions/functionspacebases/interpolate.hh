// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH

#include <memory>
#include <vector>

#include <dune/common/bitsetvector.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/common/functionfromcallable.hh>

namespace Dune {
namespace Functions {



/**
 * \brief Interpolate given function in discrete function space
 *
 * Notice that this function does currently only work for scalar
 * valued functions.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 */
template <class B, class C, class F>
void interpolate(const B& basis, C& coeff, F&& f)
{
  using GridView = typename B::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using FiniteElement = typename B::LocalView::Tree::FiniteElement;
  using FunctionBaseClass = typename Dune::LocalFiniteElementFunctionBase<FiniteElement>::type;

  using Range = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
  using LocalDomain = typename Element::Geometry::LocalCoordinate;


  auto&& gridView = basis.gridView();

  auto gf = makeGridViewFunction(std::forward<F>(f), gridView);

  auto localF = localFunction(gf);

  using FunctionFromCallable = typename Dune::Functions::FunctionFromCallable<Range(LocalDomain), decltype(localF), FunctionBaseClass>;

  auto basisIndexSet = basis.indexSet();
  coeff.resize(basisIndexSet.size());

  typename Dune::BitSetVector<1> processed(basisIndexSet.size(), false);
  std::vector<Range> interpolationValues;

  auto localView = basis.localView();
  auto localIndexSet = basisIndexSet.localIndexSet();

  for (const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);
    localF.bind(e);

    const auto& fe = localView.tree().finiteElement();

    // check if all components have already been processed
    bool allProcessed = true;
    for (size_t i=0; i<fe.localBasis().size(); ++i)
      allProcessed = allProcessed and processed[localIndexSet.index(i)[0]][0];

    if (not(allProcessed))
    {
      // We cannot use localF directly because interpolate requires a Dune::VirtualFunction like interface
      fe.localInterpolation().interpolate(FunctionFromCallable(localF), interpolationValues);
      for (size_t i=0; i<fe.localBasis().size(); ++i)
      {
        size_t index = localIndexSet.index(i)[0];
        if (not(processed[index][0]))
          coeff[index] = interpolationValues[i];
        processed[index][0] = true;
      }
    }
  }
}



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
