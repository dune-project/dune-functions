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

template <class GlobalBasis, std::size_t I>
void test(GlobalBasis const& basis, index_constant<I>)
{
  std::cout << "basis " << I << ":\n";
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

  using Tester = TaylorHoodBases<2>;
  Tester tester{};
  Hybrid::forEach(range(index_constant<Tester::num_bases>{}),
    [&tester](auto const _i) {
      test(tester.basis(_i), _i);
    });

  return 0;
}