// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/taylorhoodbasis.hh>

using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);
  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  typedef TaylorHoodBasis<GridView> Basis;

  const GridView& gridView = grid.leafGridView();
  Basis feBasis(gridView);

  typedef Basis::MultiIndex MultiIndex;

  // Sample the function f(x,y) = x on the grid vertices
  // If we use that as the coefficients of a finite element function,
  // we know its integral and can check whether quadrature returns
  // the correct result
  std::array<std::vector<double>,dim> x;
  for (unsigned int i = 0; i<dim; i++)
    x[i].resize(feBasis.size({i}));

  // TODO: Implement interpolation properly using the global basis.
  for (auto it = gridView.begin<dim>(); it != gridView.end<dim>(); ++it)
    x[1][gridView.indexSet().index(*it)] = it->geometry().corner(0)[0];

  // Objects required in the local context
  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.localIndexSet();
  std::vector<double> coefficients(localView.maxSize());

  // Loop over elements and integrate over the function
  double integral = 0;
  for (auto it = gridView.begin<0>(); it != gridView.end<0>(); ++it)
  {
    localView.bind(*it);
    localIndexSet.bind(localView);

    // paranoia checks
    assert(localView.size() == localIndexSet.size());
    assert(&(localView.globalBasis()) == &(feBasis));
    assert(&(localIndexSet.localView()) == &(localView));

    // copy data from global vector
    coefficients.resize(localIndexSet.size());
    for (size_t i=0; i<localIndexSet.size(); i++)
    {
      MultiIndex mi = localIndexSet.index(i);
      coefficients[i] = x[mi[0]][mi[1]];
    }

    // get access to the finite element
    using namespace Dune::TypeTree::Indices;
    typedef Basis::LocalView::Tree Tree;
    auto& p_leaf = TypeTree::child(localView.tree(),_1);
    auto& v_leaf = TypeTree::child(localView.tree(),_0,1);
    auto& localFiniteElement = p_leaf.finiteElement();

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
        integral += coefficients[p_leaf.localIndex(i)] * shapeFunctionValues[i] * quad[pt].weight() * integrationElement;
      }
    }

    // unbind
    localIndexSet.unbind();
    localView.unbind();
  }

  std::cout << "Computed integral is " << integral << std::endl;
  assert(std::abs(integral-0.5)< 1e-10);

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
