// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_TEST_CONTINUITYCHECKS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_TEST_CONTINUITYCHECKS_HH


using namespace Dune;

namespace Dune::Functions::Test
{
  /** \brief Check whether a given (vector-)function has continuous normal components at the element boundaries
   *
   * This is a necessary requirement for H(div)-conforming finite element functions
   *
   * \param f The function to be checked
   */
  template<class F>
  bool checkNormalContinuity(const F& f)
  {
    using GridView = typename F::Basis::GridView;
    using ctype = typename GridView::ctype;

    GridView gridView = f.basis().gridView();

    bool passed = true;

    auto insideLocalFunction = localFunction(f);
    auto outsideLocalFunction = localFunction(f);

    // Loop over all intersections
    for (const auto& element : elements(gridView))
    {
      for (const auto& intersection : intersections(gridView, element))
      {
        // Skip boundary intersections
        if (!intersection.neighbor())
          continue;

        insideLocalFunction.bind(intersection.inside());
        outsideLocalFunction.bind(intersection.outside());

        // Use quadrature rule to produce test points on the intersection
        const auto& quadRule = QuadratureRules<ctype,decltype(intersection)::mydimension>::rule(intersection.type(), 4);

        for (const auto& quadPoint : quadRule)
        {
          const auto& testPoint = quadPoint.position();

          auto insideValue  = insideLocalFunction(intersection.geometryInInside().global(testPoint));
          auto outsideValue = outsideLocalFunction(intersection.geometryInOutside().global(testPoint));

          auto normalJump = (insideValue-outsideValue) * intersection.unitOuterNormal(testPoint);

          if (std::abs(normalJump) > 1e-10)
          {
            std::cout << "Normal contribution is not continuous at " << testPoint << std::endl;
            passed = false;
          }
        }
      }
    }

    return passed;
  }

}    // namespace Dune::Functions::Test

#endif
