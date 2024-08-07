// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/taylorhoodbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

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


  test.subTest(checkBasis(feBasis, EnableContinuityCheck()));


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
  std::vector<double> coefficients(localView.maxSize());

  // Loop over elements and integrate over the function
  double integral = 0;
  for (auto it = gridView.begin<0>(); it != gridView.end<0>(); ++it)
  {
    localView.bind(*it);

    // paranoia checks
    assert(&(localView.globalBasis()) == &(feBasis));

    // copy data from global vector
    coefficients.resize(localView.size());
    for (size_t i=0; i<localView.size(); i++)
    {
      MultiIndex mi = localView.index(i);
      coefficients[i] = x[mi[0]][mi[1]];
    }

    // get access to the finite element
    using namespace Dune::Indices;
    using Tree [[maybe_unused]] = Basis::LocalView::Tree;
    auto& p_leaf = TypeTree::child(localView.tree(),_1);
    [[maybe_unused]] auto& v_leaf = TypeTree::child(localView.tree(),_0,1);
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
    localView.unbind();
  }

  std::cout << "Computed integral is " << integral << std::endl;
  assert(std::abs(integral-0.5)< 1e-10);

  /////////////////////////////////////////////////////////////////
  //  Standard tests
  /////////////////////////////////////////////////////////////////

  // check TaylorHoodBasis created 'manually'
  {
    typedef GridType::LeafGridView GridView;
    const GridView& gridView = grid.leafGridView();
    Functions::TaylorHoodBasis<GridView> basis(gridView);
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  // check RaviartThomasBasis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid.leafGridView(), taylorHood());
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  return test.exit();

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
