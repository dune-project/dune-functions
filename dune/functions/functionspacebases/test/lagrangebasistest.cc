// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

using namespace Dune;
using namespace Dune::Functions;


// Create grid with 2 cubes in 2d or 3d. While the first one is
// oriented according to the reference element, the second one
// is twisted, such that we would get a non-conforming basis if
// face DOFs are not permuted properly.
template<class Grid>
auto createNonUniformCubeGrid()
{
  Dune::GridFactory<Grid> factory;

  if constexpr (Grid::dimension == 2)
  {
    factory.insertVertex({0., 0.});
    factory.insertVertex({1., 0.});
    factory.insertVertex({0., 1.});
    factory.insertVertex({1., 1.});
    factory.insertVertex({2., 0.});
    factory.insertVertex({2., 1.});
    factory.insertElement(Dune::GeometryTypes::cube(2), {0, 1, 2, 3});
    factory.insertElement(Dune::GeometryTypes::cube(2), {5, 3, 4, 1});
  }

  if constexpr (Grid::dimension == 3)
  {
    factory.insertVertex({0., 0., 0.});
    factory.insertVertex({1., 0., 0.});
    factory.insertVertex({0., 1., 0.});
    factory.insertVertex({1., 1., 0.});
    factory.insertVertex({0., 0., 1.});
    factory.insertVertex({1., 0., 1.});
    factory.insertVertex({0., 1., 1.});
    factory.insertVertex({1., 1., 1.});
    factory.insertVertex({2., 0., 0.});
    factory.insertVertex({2., 1., 0.});
    factory.insertVertex({2., 0., 1.});
    factory.insertVertex({2., 1., 1.});
    factory.insertElement(Dune::GeometryTypes::cube(3), {0, 1, 2, 3, 4, 5, 6, 7});
    factory.insertElement(Dune::GeometryTypes::cube(3), {10, 11, 5, 7, 8, 9, 1, 3});
  }

  return factory.createGrid();
}


template<class Basis>
double benchmarkBind(const Basis& basis)
{
  auto localView = basis.localView();
  auto timer = Dune::Timer();
  for(const auto& element : elements(basis.gridView()))
    localView.bind(element);
  return timer.elapsed();
}



int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;



  using namespace Dune::Functions::BasisFactory;

  {
    const int dim = 2;
    using Grid = Dune::UGGrid<dim>;

    Dune::GridFactory<Grid> factory;
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({1,0});
    factory.insertVertex({1,1});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {0,1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {1,2,3});

    std::unique_ptr<Grid> grid(factory.createGrid());

    grid->globalRefine(2);

    auto gridView = grid->leafGridView();
    auto basis = makeBasis(gridView, lagrange<3>());

    test.subTest(checkBasis(basis, EnableContinuityCheck()));

    std::vector<double> v;
    v.resize(basis.size(), 0);
    v[5] = 1;
    auto v_f = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis,v);

    SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, Dune::refinementLevels(5));
    vtkWriter.addVertexData(v_f, VTK::FieldInfo("lambda_5", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("debug");

    // Modify grid, update basis and check again
    const auto firstEntity = gridView.template begin<0>();
    grid->mark(1, *firstEntity);
    grid->adapt();
    basis.update(grid->leafGridView());

    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  {
    auto gridPtr = createNonUniformCubeGrid<UGGrid<2>>();
    auto& grid = *gridPtr;
    grid.globalRefine(2);

    // Polynomial orders to check
    auto orders = std::index_sequence<1,2,3,4,5>{};

    // Create bases with given orders, once provided as
    // compile-time parameter and once as run-time parameter
    auto bases = Dune::unpackIntegerSequence([&](auto... order) {
      return Dune::TupleVector(
        makeBasis(grid.leafGridView(), lagrange<order>())...,
        makeBasis(grid.leafGridView(), lagrange(order))...
      );
    }, orders);

    // Check each basis
    Dune::Hybrid::forEach(bases, [&](auto& basis)
    {
      test.subTest(checkBasis(basis, EnableContinuityCheck()));
    });

    // Refine grid adaptively
    const auto firstEntity = grid.leafGridView().template begin<0>();
    grid.mark(1, *firstEntity);
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    // Update bases and check again
    Dune::Hybrid::forEach(bases, [&](auto& basis)
    {
      basis.update(grid.leafGridView());
      test.subTest(checkBasis(basis, EnableContinuityCheck()));
    });

    // Now modify the grid, update basis and check again
//     auto basis1CT = makeBasis(grid.leafGridView(), lagrange<1>());

//     grid.globalRefine(7);
//     std::cout << "2d cube, order 1 " << benchmarkBind(makeBasis(grid.leafGridView(), lagrange<1>())) << std::endl;
//     std::cout << "2d cube, order 2 " << benchmarkBind(makeBasis(grid.leafGridView(), lagrange<2>())) << std::endl;
//     std::cout << "2d cube, order 3 " << benchmarkBind(makeBasis(grid.leafGridView(), lagrange<3>())) << std::endl;
//     std::cout << "2d cube, order 4 " << benchmarkBind(makeBasis(grid.leafGridView(), lagrange<4>())) << std::endl;
//     std::cout << "2d cube, order 5 " << benchmarkBind(makeBasis(grid.leafGridView(), lagrange<5>())) << std::endl;
//     std::cout << "2d cube, order 6 " << benchmarkBind(makeBasis(grid.leafGridView(), lagrange<6>())) << std::endl;
  }

  {
    auto gridPtr = StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0., 0.}, {1., 1.}, {1, 1});
    auto& grid = *gridPtr;
    grid.globalRefine(2);

    // Polynomial orders to check
    auto orders = std::index_sequence<1,2,3,4,5>{};

    // Create bases with given orders, once provided as
    // compile-time parameter and once as run-time parameter
    auto bases = Dune::unpackIntegerSequence([&](auto... order) {
      return Dune::TupleVector(
        makeBasis(grid.leafGridView(), lagrange<order>())...,
        makeBasis(grid.leafGridView(), lagrange(order))...
      );
    }, orders);

    // Check each basis
    Dune::Hybrid::forEach(bases, [&](auto& basis)
    {
      test.subTest(checkBasis(basis, EnableContinuityCheck()));
    });

    // Refine grid adaptively
    const auto firstEntity = grid.leafGridView().template begin<0>();
    grid.mark(1, *firstEntity);
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    // Update bases and check again
    Dune::Hybrid::forEach(bases, [&](auto& basis)
    {
      basis.update(grid.leafGridView());
      test.subTest(checkBasis(basis, EnableContinuityCheck()));
    });

  }

  {
    auto gridPtr = createNonUniformCubeGrid<UGGrid<3>>();
    auto& grid = *gridPtr;
    grid.globalRefine(2);

    // Polynomial orders to check
    auto orders = std::index_sequence<1,2,3,4,5>{};

    // Create bases with given orders, once provided as
    // compile-time parameter and once as run-time parameter
    auto bases = Dune::unpackIntegerSequence([&](auto... order) {
      return Dune::TupleVector(
        makeBasis(grid.leafGridView(), lagrange<order>())...,
        makeBasis(grid.leafGridView(), lagrange(order))...
      );
    }, orders);

    // Check each basis
    Dune::Hybrid::forEach(bases, [&](auto& basis)
    {
      test.subTest(checkBasis(basis, EnableContinuityCheck()));
    });

    // Since pyramid and prism elements are only implemented
    // For orders up to 2, we only check local adaptation with these
    // compile-time parameter and once as run-time parameter
    auto adaptOrders = std::index_sequence<1,2>{};

    auto adaptBases = Dune::unpackIntegerSequence([&](auto... order) {
      return Dune::TupleVector(
        makeBasis(grid.leafGridView(), lagrange<order>())...,
        makeBasis(grid.leafGridView(), lagrange(order))...
      );
    }, adaptOrders);

    // Refine grid adaptively
    const auto firstEntity = grid.leafGridView().template begin<0>();
    grid.mark(1, *firstEntity);
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    // Update bases and check again
    Dune::Hybrid::forEach(adaptBases, [&](auto& basis)
    {
      basis.update(grid.leafGridView());
      test.subTest(checkBasis(basis, EnableContinuityCheck()));
    });

  }

  {
    auto gridPtr = StructuredGridFactory<UGGrid<3>>::createSimplexGrid({0., 0., 0.}, {1., 1., 1.}, {1, 1, 1});
    auto& grid = *gridPtr;
    grid.globalRefine(2);

    // Polynomial orders to check
    auto orders = std::index_sequence<1,2,3,4,5>{};

    // Create bases with given orders, once provided as
    // compile-time parameter and once as run-time parameter
    auto bases = Dune::unpackIntegerSequence([&](auto... order) {
      return Dune::TupleVector(
        makeBasis(grid.leafGridView(), lagrange<order>())...,
        makeBasis(grid.leafGridView(), lagrange(order))...
      );
    }, orders);

    // Check each basis
    Dune::Hybrid::forEach(bases, [&](auto& basis)
    {
      test.subTest(checkBasis(basis, EnableContinuityCheck()));
    });

    // Refine grid adaptively
    const auto firstEntity = grid.leafGridView().template begin<0>();
    grid.mark(1, *firstEntity);
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    // Update bases and check again
    Dune::Hybrid::forEach(bases, [&](auto& basis)
    {
      basis.update(grid.leafGridView());
      test.subTest(checkBasis(basis, EnableContinuityCheck()));
    });

  }

  {
    std::unique_ptr<OneDGrid> grid
      = StructuredGridFactory<OneDGrid>::createCubeGrid({0}, {1}, {10});

    auto gridView = grid->leafGridView();

    auto basis3 = makeBasis(gridView, lagrange<3>());
    auto basis2 = makeBasis(gridView, lagrange(2));
    auto basis3f = makeBasis(gridView, lagrange<3,float>());
    auto basis2f = makeBasis(gridView, lagrange<float>(2));

    test.subTest(checkBasis(basis3, EnableContinuityCheck()));
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
    test.subTest(checkBasis(basis3f, EnableContinuityCheck()));
    test.subTest(checkBasis(basis2f, EnableContinuityCheck()));

    // Modify grid, update basis and check again
    const auto firstEntity = gridView.template begin<0>();
    grid->mark(1, *firstEntity);
    grid->adapt();

    basis3.update(grid->leafGridView());
    basis2.update(grid->leafGridView());
    basis3f.update(grid->leafGridView());
    basis2f.update(grid->leafGridView());

    test.subTest(checkBasis(basis3, EnableContinuityCheck()));
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
    test.subTest(checkBasis(basis3f, EnableContinuityCheck()));
    test.subTest(checkBasis(basis2f, EnableContinuityCheck()));
  }


  return test.exit();
}
