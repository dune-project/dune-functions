// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>

//#include <dune/common/exceptions.hh>

//#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/common/localfunction.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>





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


  [[maybe_unused]] auto entitySet = EntitySet(grid.leafGridView());


  GridFunction<double(FieldVector<double,dim>), EntitySet> f;

  // Since we default-constructed the function it must be invalid.
  if (f)
  {
    std::cerr << "Default-constructed GridFunction object is not invalid!\n";
    return 1;
  }

  // Assign a test function to the GridFunction object
  auto testFunction = [](FieldVector<double,dim> x) { return x[0]*x[1]; } ;
  f = Functions::makeGridViewFunction(testFunction, grid.leafGridView());

  // Since we assigned something to f it must be valid now.
  if (!f)
  {
    std::cerr << "Nonempty GridFunction object is not valid!\n";
    return 1;
  }

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
