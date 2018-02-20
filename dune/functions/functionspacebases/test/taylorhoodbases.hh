#pragma once

#include <array>
#include <memory>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

namespace Dune { namespace Functions
{
  template <int dim, std::size_t K = 1>
  class TaylorHoodBases
  {
    using GridType = YaspGrid<dim>;

  public:
    static const std::size_t num_bases = 7;
    static const std::size_t num_false_bases = 1;

    TaylorHoodBases()
    {
      FieldVector<double,dim> bbox = {1, 1};
      std::array<int,dim> num = {1, 1};

      grid_.reset(new GridType(bbox, num));
    }

    GridType& grid() { return *grid_; }
    auto gridView() { return grid_->leafGridView(); }

    auto basis(index_constant<0>) // Root: blockedLexicographic, Velocity: flatLexicographic
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

    auto basis(index_constant<1>) // Root: blockedLexicographic, Velocity: flatInterleaved
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

    auto basis(index_constant<2>) // Root: blockedLexicographic, Velocity: blockedLexicographic
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

    auto basis(index_constant<3>) // Root: blockedLexicographic, Velocity: leafBlockedInterleaved
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

    auto basis(index_constant<4>) // Root: flatLexicographic, Velocity/Pressure: flatLexicographic
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

    auto basis(index_constant<5>) // Root: blockedLexicographic, Velocity/Pressure: blockedLexicographic
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
          blockedLexicographic()
        ));
    }

    auto basis(index_constant<6>) // Root: blockedLexicographic, Velocity/Pressure: on the same level
    {
      using namespace Dune::Functions::BasisBuilder;
      return makeBasis(
        gridView(),
        composite(
          lagrange<K+1>(),
          lagrange<K+1>(),
          lagrange<K>(),
          blockedLexicographic()
        ));
    }

    auto false_basis(index_constant<0>) // Root: flatLexicographic, Velocity/Pressure: blockedLexicographic
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

  private:
    std::unique_ptr<GridType> grid_;
  };

}} // end namespace Dune::Functions
