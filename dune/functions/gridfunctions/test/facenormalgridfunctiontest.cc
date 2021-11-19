// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stringutility.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/functions/gridfunctions/facenormalgridfunction.hh>

template<class GridView>
Dune::TestSuite checkFaceNormalGridFunction(const GridView& gridView, std::string name)
{
  using Dune::Functions::Concept::isGridFunction;
  using F = Dune::Functions::FaceNormalGridFunction<GridView>;
  using EntitySet = typename F::EntitySet;
  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = Domain;

  auto testSuite = Dune::TestSuite(name);

  // Check is GridFunction concept is satisfied
  testSuite.check(isGridFunction<F, Range(Domain), EntitySet>())
    << "FaceNormalGridFunction does not model GridFunction concept";

  // Check if normals coincide with intersection normals
  // up to given tolerance at quadrature point.
  double TOL = 1e-14;
  std::size_t sampleQuadratureOrder = 10;
  constexpr int dimension = GridView::dimension;
  auto normals = Dune::Functions::FaceNormalGridFunction(gridView);
  auto localNormals = localFunction(normals);
  for(auto&& element : elements(gridView))
  {
    localNormals.bind(element);
    for(const auto& intersection: intersections(gridView, element))
    {
      const auto& quad = Dune::QuadratureRules<double, dimension-1>::rule(intersection.type(), sampleQuadratureOrder);
      double mismatch = 0.0;
      for (const auto& quadPoint : quad)
      {
        auto normal = intersection.unitOuterNormal(quadPoint.position());
        normal -= localNormals(intersection.geometryInInside().global(quadPoint.position()));
        mismatch = std::max(mismatch, normal.infinity_norm());
      }
      testSuite.check(mismatch <= TOL)
        << "Normal mismatch of " << Dune::formatString("%12.5e", mismatch)
        << " for face " << intersection.indexInInside()
        << " of " << intersection.inside().type()
        << " with index " << gridView.indexSet().index(intersection.inside());
    }
  }

  return testSuite;
}



int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite;

  {
    using Grid = Dune::YaspGrid<2>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid({{0,0}}, {{1,1}}, {{2,2}});
    grid->globalRefine(2);
    testSuite.subTest(checkFaceNormalGridFunction(grid->leafGridView(), "YaspGrid<2>"));
  }

  {
    using Grid = Dune::YaspGrid<3>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid({{0,0,0}}, {{1,1,1}}, {{2,2,2}});
    grid->globalRefine(2);
    testSuite.subTest(checkFaceNormalGridFunction(grid->leafGridView(), "YaspGrid<3>"));
  }

  {
    using Grid = Dune::UGGrid<2>;
    auto grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid({{0,0}}, {{1,1}}, {{1,1}});
    grid->globalRefine(2);
    testSuite.subTest(checkFaceNormalGridFunction(grid->leafGridView(), "UGGrid<2> (triangles)"));
  }

  {
    using Grid = Dune::UGGrid<2>;
    auto factory = Dune::GridFactory<Grid>();
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({1,0});
    factory.insertVertex({2,2});
    factory.insertElement(Dune::GeometryTypes::cube(2), {0,1,2,3});
    auto grid = factory.createGrid();
    grid->globalRefine(2);
    testSuite.subTest(checkFaceNormalGridFunction(grid->leafGridView(), "UGGrid<2> (nonaffine rectangles)"));
  }

  {
    using Grid = Dune::UGGrid<3>;
    auto grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid({{0,0,0}}, {{1,1,1}}, {{1,1,1}});
    grid->globalRefine(2);
    testSuite.subTest(checkFaceNormalGridFunction(grid->leafGridView(), "UGGrid<2> (triangles)"));
  }

  {
    using Grid = Dune::UGGrid<3>;
    auto factory = Dune::GridFactory<Grid>();
    factory.insertVertex({0,0,0});
    factory.insertVertex({1,0,0});
    factory.insertVertex({0,1,0});
    factory.insertVertex({1,1,0});
    factory.insertVertex({0,0,1});
    factory.insertVertex({1,0,1});
    factory.insertVertex({0,1,1});
    factory.insertVertex({2,2,2});
    factory.insertElement(Dune::GeometryTypes::cube(3), {0,1,2,3,4,5,6,7});
    auto grid = factory.createGrid();
    grid->globalRefine(2);
    testSuite.subTest(checkFaceNormalGridFunction(grid->leafGridView(), "UGGrid<3> (nonaffine cubes)"));
  }

  return testSuite.exit();
}
catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
