#include "config.h"

#include <dune/common/std/type_traits.hh>
#include <dune/functions/functionspacebases/matrixvectorgenerator.hh>
#include <dune/functions/functionspacebases/hierarchicresizer.hh>
#include <dune/functions/functionspacebases/sizeinfo.hh>

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

template <class BlockingTag>
struct ContainsLeafBlocked
    : std::false_type {};

template <std::size_t N>
struct ContainsLeafBlocked<Blocking::LeafBlocked<N>>
    : std::true_type {};

template <class... B>
struct ContainsLeafBlocked<Blocking::Blocked<B...>>
    : Std::disjunction<ContainsLeafBlocked<B>...> {};


// test vector resize
template <class GlobalBasis, std::size_t I>
void test1(GlobalBasis const& basis, index_constant<I>)
{
  std::cout << "basis " << I << ":\n";
  using Blocking = Blocking_t<GlobalBasis>;

  using Vector = VectorGenerator_t<FieldVector<double,1>, BlockVector, Blocking>;
  using VectorResizer = HierarchicVectorResizer<Vector, Blocking>;
  Vector v;
  VectorResizer vectorResizer(v);
  vectorResizer(sizeInfo(basis));
  std::cout << "v.size = " << hybridSize(v) << "\n";
}


template <class RowBasis, class ColBasis, std::size_t I, std::size_t J>
void test2_impl(RowBasis const& rowBasis, ColBasis const& colBasis, index_constant<I>, index_constant<J>, std::true_type)
{
  /* Mixed basis for rows and columns not allowed, if one is LeafBlocked */
}

// test matrix resize
template <class RowBasis, class ColBasis, std::size_t I, std::size_t J>
void test2_impl(RowBasis const& rowBasis, ColBasis const& colBasis, index_constant<I>, index_constant<J>, std::false_type)
{
  std::cout << "rowbasis " << I << ", colbasis " << J << ":\n";
  using RowBlocking = Blocking_t<RowBasis>;
  using ColBlocking = Blocking_t<ColBasis>;

  using Matrix = MatrixGenerator_t<FieldMatrix<double,1,1>, BCRSMatrix, RowBlocking, ColBlocking>;
  using MatrixResizer = HierarchicMatrixResizer<Matrix, RowBlocking, ColBlocking>;
  Matrix m;
  MatrixResizer matrixResizer(m);
  matrixResizer(sizeInfo(rowBasis), sizeInfo(colBasis));
  std::cout << "m.size = " << hybridNumRows(m) << ", " << hybridNumCols(m) << "\n";
}

template <class RowBasis, class ColBasis, std::size_t I, std::size_t J>
void test2(RowBasis const& rowBasis, ColBasis const& colBasis, index_constant<I> _i, index_constant<J> _j)
{
  test2_impl(rowBasis, colBasis, _i, _j, Std::disjunction<ContainsLeafBlocked<Blocking_t<RowBasis>>, ContainsLeafBlocked<Blocking_t<ColBasis>>>{});
}

template <class RowBasis, std::size_t I>
void test2(RowBasis const& rowBasis, RowBasis const&, index_constant<I> _i, index_constant<I>)
{
  test2_impl(rowBasis, rowBasis, _i, _i, std::false_type{});
}


int main(int argc, char** argv)
{
  MPIHelper::instance(argc, argv);

  std::cout << "test1...\n";
  using Tester = TaylorHoodBases<2>;
  Tester tester{};
  Hybrid::forEach(range(index_constant<Tester::num_bases>{}),
    [&tester](auto const _i) {
      test1(tester.basis(_i), _i);
    });

  std::cout << "test2...\n";
  Hybrid::forEach(range(index_constant<Tester::num_bases>{}),
    [&](auto const _i) {
      Hybrid::forEach(range(index_constant<Tester::num_bases>{}),
        [&](auto const _j) {
          test2(tester.basis(_i), tester.basis(_j), _i, _j);
        });
    });
  return 0;
}