// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

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

    test.subTest(checkBasis(basis));

    std::vector<double> v;
    v.resize(basis.size(), 0);
    v[5] = 1;
    auto v_f = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis,v);

    SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, Dune::refinementLevels(5));
    vtkWriter.addVertexData(v_f, VTK::FieldInfo("lambda_5", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("debug");

  }

  {
    std::unique_ptr<OneDGrid> grid
      = StructuredGridFactory<OneDGrid>::createCubeGrid({0}, {1}, {10});

    auto gridView = grid->levelGridView(0);

    {
      auto basis = makeBasis(gridView, lagrange<3>());
      test.subTest(checkBasis(basis));
    }

    {
      auto basis = makeBasis(gridView, lagrange(2));
      test.subTest(checkBasis(basis));
    }

    {
      auto basis = makeBasis(gridView, lagrange<3,float>());
      test.subTest(checkBasis(basis));
    }

    {
      auto basis = makeBasis(gridView, lagrange<float>(2));
      test.subTest(checkBasis(basis));
    }

  }


  return test.exit();
}
