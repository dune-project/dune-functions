#include "config.h"

#include <dune/functions/functionspacebases/matrixvectorgenerator.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#else
#error "Need dune-istl for this test"
#endif

#include "taylorhoodbases.hh"

using namespace Dune;
using namespace Dune::Functions;

template <class T>
void printType()
{
  std::cout << __PRETTY_FUNCTION__ << "\n";
}

template <class T>
void printType(T const&)
{
  std::cout << __PRETTY_FUNCTION__ << "\n";
}


void test1(TaylorHoodBases& tester)
{
  auto basis = tester.basis1();
  using GlobalBasis = decltype(basis);
  using Blocking = Blocking_t<GlobalBasis>;
  printType<Blocking>();

  using Vector = VectorGenerator_t<FieldVector<double,1>, BlockVector, Blocking>;
  printType<Vector>();

  using Matrix = MatrixGenerator_t<FieldMatrix<double,1,1>, BCRSMatrix, Blocking, Blocking>;
  printType<Matrix>();
}

void test2(TaylorHoodBases& tester)
{
  auto basis = tester.basis2();
  using GlobalBasis = decltype(basis);
  using Blocking = Blocking_t<GlobalBasis>;
  printType<Blocking>();

  using Vector = VectorGenerator_t<FieldVector<double,1>, BlockVector, Blocking>;
  printType<Vector>();

  using Matrix = MatrixGenerator_t<FieldMatrix<double,1,1>, BCRSMatrix, Blocking, Blocking>;
  printType<Matrix>();
}

void test3(TaylorHoodBases& tester)
{
  auto basis = tester.basis3();
  using GlobalBasis = decltype(basis);
  using Blocking = Blocking_t<GlobalBasis>;
  printType<Blocking>();

  using Vector = VectorGenerator_t<FieldVector<double,1>, BlockVector, Blocking>;
  printType<Vector>();

  using Matrix = MatrixGenerator_t<FieldMatrix<double,1,1>, BCRSMatrix, Blocking, Blocking>;
  printType<Matrix>();
}

void test4(TaylorHoodBases& tester)
{
  auto basis = tester.basis4();
  using GlobalBasis = decltype(basis);
  using Blocking = Blocking_t<GlobalBasis>;
  printType<Blocking>();

  using Vector = VectorGenerator_t<FieldVector<double,1>, BlockVector, Blocking>;
  printType<Vector>();

  using Matrix = MatrixGenerator_t<FieldMatrix<double,1,1>, BCRSMatrix, Blocking, Blocking>;
  printType<Matrix>();
}

void test5(TaylorHoodBases& tester)
{}

void test6(TaylorHoodBases& tester)
{
  auto basis = tester.basis6();
  using GlobalBasis = decltype(basis);
  using Blocking = Blocking_t<GlobalBasis>;
  printType<Blocking>();

  using Vector = VectorGenerator_t<FieldVector<double,1>, BlockVector, Blocking>;
  printType<Vector>();

  using Matrix = MatrixGenerator_t<FieldMatrix<double,1,1>, BCRSMatrix, Blocking, Blocking>;
  printType<Matrix>();
}

int main(int argc, char** argv)
{
  MPIHelper::instance(argc, argv);

  TaylorHoodBases tester{};
  test1(tester);
  test2(tester);
  test3(tester);
  test4(tester);
  test5(tester);
  test6(tester);

  return 0;
}