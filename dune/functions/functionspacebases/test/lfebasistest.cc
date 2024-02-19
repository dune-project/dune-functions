// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/lfeprebasismixin.hh>
#include <dune/functions/functionspacebases/refinedlagrangebasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;
using namespace Dune::Functions;

// define a pre-basis based on the LFEPreBasisMixin
template <class GV, class R = double>
class RefinedP0PreBasis :
  public Dune::Functions::LFEPreBasisMixin<GV, Dune::RefinedP0LocalFiniteElement<typename GV::ctype,R,GV::dimension>>
{
  using LFE = Dune::RefinedP0LocalFiniteElement<typename GV::ctype,R,GV::dimension>;
  using Base = LFEPreBasisMixin<GV, LFE>;
  static const int dim = GV::dimension;
public:
  RefinedP0PreBasis (const GV& gv) :
    Base(gv, [](Dune::GeometryType gt, int) { return (gt.dim() == dim) ? (1 << dim) : 0; })
  {}
};


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  using namespace Dune::Functions::BasisFactory;

  {
    const int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    using Factory = Dune::StructuredGridFactory<Grid>;

    auto gridPtr = Factory::createSimplexGrid({0.0,0.0},{1.0,1.0},{1u,1u});
    gridPtr->globalRefine(2);

    auto gridView = gridPtr->leafGridView();
    using GridView = decltype(gridView);

    // compare two bases that should be exactly identical, but that are implemented differently
    auto basis0 = makeBasis(gridView, refinedLagrange<0>());
    auto basis1 = DefaultGlobalBasis<RefinedP0PreBasis<GridView>>(gridView);

    test.subTest(checkBasis(basis0));
    test.subTest(checkBasis(basis1));

    auto localView0 = basis0.localView();
    auto localView1 = basis1.localView();

    test.check(localView0.maxSize() == localView1.maxSize(), "maxSize");

    for (auto const& e : elements(gridView))
    {
      localView0.bind(e);
      localView1.bind(e);

      test.check(localView0.size() == localView1.size(), "size");
      for (std::size_t i = 0; i < localView0.size(); ++i)
      {
        test.check(localView0.index(i) == localView1.index(i), "index");
      }

      auto const& lfe0 = localView0.tree().finiteElement();
      auto const& lfe1 = localView1.tree().finiteElement();

      using LFE0 = std::decay_t<decltype(lfe0)>;
      using LFE1 = std::decay_t<decltype(lfe1)>;
      static_assert(std::is_same_v<LFE0,LFE1>);
    }
  }

  return test.exit();
}
