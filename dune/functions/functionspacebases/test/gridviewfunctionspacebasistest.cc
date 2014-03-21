// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/pq1nodalbasis.hh>

using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[]) try
{
  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
//  FieldVector<double,dim> l = {{21.0, 4.0}};
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);

  // Test whether PQ1FunctionSpaceBasis.hh can be instantiated on the leaf view
  typedef GridType::LeafGridView GridView;
  typedef PQ1NodalBasis<GridView> Basis;

  const GridView& gridView = grid.leafGridView();
  Basis feBasis(gridView);

  typedef Basis::MultiIndex MultiIndex;

  // Sample the function f(x,y) = x on the grid vertices
  // If we use that as the coefficients of a finite element function,
  // we know its integral and can check whether quadrature returns
  // the correct result
  std::vector<double> x(feBasis.subIndexCount());



  // TODO: Implement interpolation properly using the global basis.
  for (auto it = gridView.begin<dim>(); it != gridView.end<dim>(); ++it)
    x[gridView.indexSet().index(*it)] = it->geometry().corner(0)[0];

  // Loop over elements and integrate over the function
  double integral = 0;
  std::vector<MultiIndex> globalIndices(feBasis.maxLocalSize());
  for (auto it = gridView.begin<0>(); it != gridView.end<0>(); ++it)
  {
    auto localView = feBasis.localView();
    localView.bind(*it);

    typedef Basis::LocalView::Tree Tree;
    auto& treeImp = localView.tree();
    const Tree::Interface& tree = treeImp;

    treeImp.generateMultiIndices(globalIndices.begin());

    auto& localFiniteElement = tree.finiteElement();

    // A quadrature rule
    const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(it->type(), 1);

    // Loop over all quadrature points
    for ( size_t pt=0; pt < quad.size(); pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = it->geometry().integrationElement(quadPos);

      // Evaluate all shape function values at this point
      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

      // Actually compute the vector entries
      for (size_t i=0; i<localFiniteElement.localBasis().size(); i++)
      {
        assert(globalIndices[i]==tree.globalIndex(i));
        integral += x[globalIndices[i][0]] * shapeFunctionValues[i] * quad[pt].weight() * integrationElement;
      }
    }
  }

  assert(std::abs(integral-0.5)< 1e-10);
  std::cout << "Computed integral is " << integral << std::endl;

  return 0;

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
}
