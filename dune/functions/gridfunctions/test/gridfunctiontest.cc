// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

//#include <dune/common/exceptions.hh>

//#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/common/localfunction.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

#include <dune/functions/analyticfunctions/polynomial.hh>
#include <dune/functions/analyticfunctions/trigonometricfunction.hh>





using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
//  FieldVector<double,dim> l = {{21.0, 4.0}};
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);

  // Test whether PQ1FunctionSpaceBasis.hh can be instantiated on the leaf view
  using GridView = GridType::LeafGridView;
  using EntitySet = GridViewEntitySet<GridView, 0>;


  auto entitySet = EntitySet(grid.leafGridView());


  auto p = Polynomial<double>({0, 1, 2});

  GridFunction<double(double), EntitySet> f;





  return 0;

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
