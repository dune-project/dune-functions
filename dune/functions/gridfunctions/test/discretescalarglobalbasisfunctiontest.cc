// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

// #include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/functions/functionspacebases/pq2nodalbasis.hh>
#include <dune/functions/gridfunctions/new_discretescalarglobalbasisfunction.hh>


#include <dune/functions/gridfunctions/new_gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>


using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[]) try
{
  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);

  // Test whether PQ1FunctionSpaceBasis.hh can be instantiated on the leaf view
  typedef GridType::LeafGridView GridView;
//  typedef PQ1NodalBasis<GridView> Basis;
  typedef PQ2NodalBasis<GridView> Basis;

  const GridView& gridView = grid.leafGridView();
  Basis feBasis(gridView);
  auto indexSet = feBasis.indexSet();

  // Sample the function f(x,y) = x on the grid vertices
  // If we use that as the coefficients of a finite element function,
  // we know its integral and can check whether quadrature returns
  // the correct result
  std::vector<double> x(indexSet.size());

  // TODO: Implement interpolation properly using the global basis.
  for (auto it = gridView.begin<dim>(); it != gridView.end<dim>(); ++it)
    x[gridView.indexSet().index(*it)] = it->geometry().corner(0)[0];

  // generate a discrete function to evaluate the integral
  Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasis),decltype(x)> fOriginal(feBasis,x);

  using EntitySet = decltype(fOriginal)::EntitySet;
  using GridView = decltype(fOriginal)::GridView;
  using Range = decltype(fOriginal)::Range;
  using Domain = decltype(fOriginal)::Domain;


//  auto& f = fOriginal;
//  GridFunction<Range(Domain), EntitySet ff = fOriginal;
  GridFunction<Range(Domain), EntitySet> f = fOriginal;




//  auto localFunction = Dune::Functions::localFunction(f);
  auto fLocal = localFunction(f);

  // Loop over elements and integrate over the function
  double integral = 0;
  for (auto it = gridView.begin<0>(); it != gridView.end<0>(); ++it)
  {
//    fLocal->bind(*it);
    fLocal.bind(*it);

    // A quadrature rule
    const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(it->type(), 1);

    // Loop over all quadrature points
    for ( size_t pt=0; pt < quad.size(); pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = it->geometry().integrationElement(quadPos);

      Dune::FieldVector<typename decltype(x)::value_type, 1> v;
      v = fLocal(quadPos);
//      fLocal->evaluate(quadPos,v);
      integral += v * quad[pt].weight() * integrationElement;
    }

//    fLocal->unbind();
    fLocal.unbind();
  }

  std::cout << "Computed integral is " << integral << std::endl;
  assert(std::abs(integral-0.5)< 1e-10);

  return 0;

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
}
