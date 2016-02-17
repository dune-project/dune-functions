// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>
#include <cmath>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/bcrsmatrix.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/raviartthomascubebasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

using namespace Dune;

//////////////////////////////////////////////////////////
//
// Test single RT0 basis functions by plotting.
// Use BasisBuilder.
//
//////////////////////////////////////////////////////////

int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {3, 3};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  using namespace Functions::BasisTags;
  using namespace Functions::BasisBuilder;

  auto basis = makeBasis(
    gridView,
    composite(
      rtcube<0>(),
      pq<0>()
    ));

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  typedef BlockVector<BlockVector<FieldVector<double,1> > > VectorType;
  typedef Dune::Functions::HierarchicVectorWrapper<VectorType, double> HierarchicVectorView;

  VectorType u;
  HierarchicVectorView(u).resize(basis);

//  using namespace TypeTree::Indices;

  for (size_t i=0; i<u[0].size(); i++) {

    u = 0;
    u[0][i][0] = 1.0;

    ////////////////////////////////////////////////////////////////////////////
    //  Make a discrete function from the FE basis and the coefficient vector
    ////////////////////////////////////////////////////////////////////////////

    using FluxRange = FieldVector<double,dim>;

    auto fluxFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<FluxRange>(basis, TypeTree::hybridTreePath(TypeTree::Indices::_0), HierarchicVectorView(u));

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Write result to VTK file
    //////////////////////////////////////////////////////////////////////////////////////////////
    SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
    vtkWriter.addVertexData(fluxFunction, VTK::FieldInfo("flux", VTK::FieldInfo::Type::vector, dim));
    vtkWriter.write("rt0-basistest-result_"+std::to_string(i+1));
  }
 }
// Error handling
 catch (Exception e) {
    std::cout << e << std::endl;
 }
