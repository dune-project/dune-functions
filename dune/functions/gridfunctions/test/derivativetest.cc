#include <config.h>
#include <dune/common/parallel/mpihelper.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>

#include <dune/istl/bvector.hh>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

#include "dune/functions/functionspacebases/interpolate.hh"
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <math.h>

#define VECTORVALUED

using namespace Dune;

// Test case:
// Define vector valued finite element function (piecewise bilinear) on three dimensional grid
// by interpolating a given nonlinear, differentiable function. Compare the derivatives of the
// approximation with the exact derivative. Repeat the same test on finer grids and check
// if the interpolation error behaves as expected.

// Compute L2 norm of the difference between a discrete global bassis function and another gridfunctions.
template <class GridView, class FunctionA, class FunctionB>
double l2difference (const GridView& gridView, FunctionA& functionA, FunctionB& functionB)
{
  double difference = 0.0;

  for (const auto& element: elements(gridView))
  {
    auto geometry = element.geometry();
    functionA.bind(element);
    functionB.bind(element);
    const auto& quad = QuadratureRules<double, GridView::dimension>::rule(element.type(), 2);

    for (const auto& quadPoint : quad)
    {
      auto quadPos = quadPoint.position();
      auto diff  = functionA(quadPos);
      diff -= functionB(quadPos);
#ifdef VECTORVALUED
      auto norm = diff.two_norm2();
#else
      auto norm = diff * diff;
#endif
      auto integrationElement = geometry.integrationElement(quadPos);
      difference += norm * quadPoint.weight() * integrationElement;
    }
  }

  difference = sqrt(difference);
  return difference;
}

// Compute L2 norm of the difference between the gradient of a discrete
// global bassis function and another gridfunctions.
template <class GridView, class FunctionA, class FunctionB>
double h1difference (const GridView& gridView, FunctionA& functionA, FunctionB& functionB)
{
  double difference = 0.0;

  for (const auto& element: elements(gridView))
  {
    auto geometry = element.geometry();
    functionA.bind(element);
    functionB.bind(element);
    const auto& quad = QuadratureRules<double, GridView::dimension>::rule(element.type(), 2);

    for (const auto& quadPoint : quad)
    {
      auto quadPos = quadPoint.position();
      auto diff = functionB(quadPos);
      // Compute the derivative of a discrete global basis function.
      auto grad = functionA.derivative(quadPos);
      for (size_t i=0; i<diff.N(); ++i)
#ifdef VECTORVALUED
        for (size_t j=0; j<diff.M(); ++j)
          diff[i][j] -= grad[i][j];
      auto norm = diff.frobenius_norm2();
#else
        diff[i] -= grad[i];
      auto norm = diff.two_norm2();
#endif
      auto integrationElement = geometry.integrationElement(quadPos);
      difference += norm * quadPoint.weight() * integrationElement;
    }
  }

  difference = sqrt(difference);
  return difference;
}

int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Define a continuous function and its derivative
  //  f(x,y,z) = [ sin(2πx) * sin(2πy) * sin(2πz) ]
  //             [ x^2 * exp(y) - z               ];
  //
  //  Df(x,y,z) = [ 2π * cos(2πx) * sin(2πy) * sin(2πz) , 2π * sin(2πx) * cos(2πy) * sin(2πz), 2π * sin(2πx) * sin(2πy) * cos(2πz) ]
  //              [ 2x * exp(y)                         , x^2 * exp(y)                       , -1                                  ];
  //
  //////////////////////////////////////////////////////////////////////////////////////////////
  auto pi = std::acos(-1.0);
  auto function = [&pi] (const auto& x)
  {
#ifdef VECTORVALUED
    FieldVector<double, 2> result(0.0);
#else
    double result(0.0);
#endif

#ifdef VECTORVALUED
    result[0] = sin(2.*pi*x[0]) * sin(2.*pi*x[1]) * sin(2.*pi*x[2]);
    result[1] = pow(x[0],2) * exp(x[1]) - x[2];
#else
    result = sin(2.*pi*x[0]) * sin(2.*pi*x[1]) * sin(2.*pi*x[2]);
#endif

    return result;
  };

  auto derivative = [&pi] (const auto& x)
  {
#ifdef VECTORVALUED
    FieldMatrix<double, 2, 3> result(0.0);

    result[0][0] = 2.*pi * cos(2.*pi*x[0]) * sin(2.*pi*x[1]) * sin(2.*pi*x[2]);
    result[0][1] = 2.*pi * sin(2.*pi*x[0]) * cos(2.*pi*x[1]) * sin(2.*pi*x[2]);
    result[0][2] = 2.*pi * sin(2.*pi*x[0]) * sin(2.*pi*x[1]) * cos(2.*pi*x[2]);

    result[1][0] = 2. * x[0] * exp(x[1]);
    result[1][1] = pow(x[0],2) * exp(x[1]);
    result[1][2] = - 1.;

#else
    FieldVector<double, 3> result(0.0);

    result[0] = 2.*pi * cos(2.*pi*x[0]) * sin(2.*pi*x[1]) * sin(2.*pi*x[2]);
    result[1] = 2.*pi * sin(2.*pi*x[0]) * cos(2.*pi*x[1]) * sin(2.*pi*x[2]);
    result[2] = 2.*pi * sin(2.*pi*x[0]) * sin(2.*pi*x[1]) * cos(2.*pi*x[2]);

#endif

    return result;
  };

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Generate a three-dimensional grid
  //////////////////////////////////////////////////////////////////////////////////////////////
  const int dim = 3;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> ur(1.0);
  std::array<int,dim> numElements = {10,10,10};
  GridType grid(ur,numElements);
  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  for (size_t refineLevel = 0; refineLevel < 3; ++refineLevel)
  {
    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Define the continuous function and its derivative as grid functions
    //////////////////////////////////////////////////////////////////////////////////////////////
    auto gridFunction = localFunction(Functions::makeGridViewFunction(function, gridView));
    auto gridDerivative = localFunction(Functions::makeGridViewFunction(derivative, gridView));

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Choose a finite element space (piecewise bilinear)
    //////////////////////////////////////////////////////////////////////////////////////////////
    using namespace Functions::BasisBuilder;
#ifdef VECTORVALUED
    auto basis = makeBasis(gridView, power<2>(pq<1>(),flatInterleaved()));
#else
    auto basis = makeBasis(gridView, pq<1>());
#endif

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Define array type formats
    //////////////////////////////////////////////////////////////////////////////////////////////
    typedef BlockVector<FieldVector<double,1> > VectorType;
    typedef Dune::Functions::HierarchicVectorWrapper<VectorType, double> HierarchicVectorView;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Define a coefficient vector corresponding to a discrete function
    //////////////////////////////////////////////////////////////////////////////////////////////
    VectorType x;
    HierarchicVectorView(x).resize(basis);
    x = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Make a discrete function from the FE basis and the coefficient vector
    //////////////////////////////////////////////////////////////////////////////////////////////
#ifdef VECTORVALUED
    using FunctionRange = FieldVector<double,2>;
#else
    using FunctionRange = double;
#endif
    auto discreteFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<FunctionRange>(basis, TypeTree::hybridTreePath(), HierarchicVectorView(x));

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Transform discrete functions to local functions
    //////////////////////////////////////////////////////////////////////////////////////////////
    auto localDiscreteFunction = localFunction(discreteFunction);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Interpolate the function
    //////////////////////////////////////////////////////////////////////////////////////////////
    interpolate(basis, Dune::TypeTree::hybridTreePath(), HierarchicVectorView(x), function);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Define VTK writer and output
    //////////////////////////////////////////////////////////////////////////////////////////////
//    SubsamplingVTKWriter<GridView> vtkWriter(gridView,1);
//#ifdef VECTORVALUED
//    vtkWriter.addVertexData(discreteFunction, VTK::FieldInfo("function", VTK::FieldInfo::Type::vector, 2));
//#else
//    vtkWriter.addVertexData(discreteFunction, VTK::FieldInfo("function", VTK::FieldInfo::Type::scalar, 1));
//#endif
//    vtkWriter.write("derivativetest");

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Check l2 difference for function and its derivative
    //////////////////////////////////////////////////////////////////////////////////////////////
    auto l2diff = l2difference(gridView, localDiscreteFunction, gridFunction);
    auto h1diff = h1difference(gridView, localDiscreteFunction, gridDerivative);

    std::cout << "Refine level: " << refineLevel << std::endl;
    std::cout << "L2 interpolation error: " << l2diff << std::endl;
    std::cout << "H1-seminorm interpolation error: " << h1diff << std::endl << std::endl;

    grid.globalRefine(1);
  }

 }
// Error handling
 catch (Exception e) {
    std::cout << e << std::endl;
 }
