#pragma once

#include <array>
#include <memory>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

namespace Dune { namespace Functions
{
  class TaylorHoodBases
  {
    static const int dim = 2;
    static const std::size_t K = 1; // pressure order for Taylor-Hood
    using GridType = YaspGrid<dim>;

  public:
    TaylorHoodBases()
    {
      FieldVector<double,dim> bbox = {1, 1};
      std::array<int,dim> num = {1, 1};

      grid_.reset(new GridType(bbox, num));
    }

    GridType& grid() { return *grid_; }
    auto gridView() { return grid_->leafGridView(); }

    auto basis1() // Root: blockedLexicographic, Velocity: flatLexicographic
    {
      using namespace Dune::Functions::BasisBuilder;
      return makeBasis(
        gridView(),
        composite(
          power<dim>(
            lagrange<K+1>(),
            flatLexicographic()),
          lagrange<K>(),
          blockedLexicographic()
        ));
    }

    auto basis2() // Root: blockedLexicographic, Velocity: flatInterleaved
    {
      using namespace Dune::Functions::BasisBuilder;
      return makeBasis(
        gridView(),
        composite(
          power<dim>(
            lagrange<K+1>(),
            flatInterleaved()),
          lagrange<K>(),
          blockedLexicographic()
        ));
    }

    auto basis3() // Root: blockedLexicographic, Velocity: blockedLexicographic
    {
      using namespace Dune::Functions::BasisBuilder;
      return makeBasis(
        gridView(),
        composite(
          power<dim>(
            lagrange<K+1>(),
            blockedLexicographic()),
          lagrange<K>(),
          blockedLexicographic()
        ));
    }

    auto basis4() // Root: blockedLexicographic, Velocity: leafBlockedInterleaved
    {
      using namespace Dune::Functions::BasisBuilder;
      return makeBasis(
        gridView(),
        composite(
          power<dim>(
            lagrange<K+1>(),
            leafBlockedInterleaved()),
          lagrange<K>(),
          blockedLexicographic()
        ));
    }

    auto basis5() // Root: flatLexicographic, Velocity/Pressure: blockedLexicographic
    {
      using namespace Dune::Functions::BasisBuilder;
      return makeBasis(
        gridView(),
        composite(
          power<dim>(
            lagrange<K+1>(),
            blockedLexicographic()),
          power<1>(
            lagrange<K>(),
            blockedLexicographic()),
          flatLexicographic()
        ));
    }

    auto basis6() // Root: flatLexicographic, Velocity/Pressure: flatLexicographic
    {
      using namespace Dune::Functions::BasisBuilder;
      return makeBasis(
        gridView(),
        composite(
          power<dim>(
            lagrange<K+1>(),
            flatLexicographic()),
          power<1>(
            lagrange<K>(),
            flatLexicographic()),
          flatLexicographic()
        ));
    }

  private:
    std::unique_ptr<GridType> grid_;
  };

}} // end namespace Dune::Functions
