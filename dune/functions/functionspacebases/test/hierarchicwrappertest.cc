// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>
#include <array>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/classname.hh>
#include <dune/common/fvector.hh>
#include <dune/common/tuplevector.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dune/typetree/utility.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/hierarchicwrapper.hh>


/**
 * \brief A Dummy size provider
 *
 * This is a mock class providing non-uniform size information.
 * It's non-uniform in the sense, that not all multi-indices are
 * do not always have the same size.
 */
template<std::size_t dim>
class HybridSizeInfoDummy
{
public:
    using size_type = std::size_t;
    using SizePrefix = Dune::ReservedVector<std::size_t, 3>;

    /**
     * \brief Construct from basis
     */
    HybridSizeInfoDummy()
    {}

    /**
     * \brief Return number possible values for next position in multi index
     */
    size_type operator()(const SizePrefix& prefix) const
    {
      return size(prefix);
    }

    /**
     * \brief Return number possible values for next position in multi index
     *
     * This shall vanish. It's just here such that this can be used
     * as size provider n place of the basis.
     */
    size_type size(const SizePrefix& prefix) const
    {
      if (prefix.size() == 0)
        return 2;
      if (prefix.size() == 1)
      {
        if (prefix[0] == 0)
          return 23;
        if (prefix[0] == 1)
          return 42;
      }
      if (prefix.size() == 2)
      {
        if (prefix[0] == 0)
          return dim;
        if (prefix[0] == 1)
          return 0;
      }
      if (prefix.size() == 3)
        return 0;
      assert(false);
    }

    operator size_type () const
    {
        return 23*dim+42;
    }

};


#if 0
template<class Vector, class SizeInfo, class SizePrefix,
  typename std::enable_if< not HasStaticSize<Vector>::value, int>::type = 0>
bool checkHierarchicVectorSize(const Vector& v, const SizeInfo& sizeInfo, SizePrefix prefix)
{
  TestSuite test;;

  test.require(v.size() == sizeInfo(SizePrefix{}))
  prefix.push_back(0);
  for (std::size_t i=0; i< v.size(); ++i)
  {
    prefix.back() = i;
    test.check(checkHierarchicVectorSize(v[i], sizeInfo, prefix)) << "Size check for entry with prefix " << prefix << " failed";
  }
  return test;
}


template<class Vector, class SizeInfo, class SizePrefix>
bool checkHierarchicVectorSize(const Vector& v, const SizeInfo& sizeInfo, SizePrefix prefix)
{
  TestSuite test;;

  test.require(v.size() == sizeInfo(SizePrefix{}))
  prefix.push_back(0);
  for (std::size_t i=0; i< v.size(); ++i)
  {
    prefix.back() = i;
    test.check(checkHierarchicVectorSize(v[i], sizeInfo, prefix)) << "Size check for entry with prefix " << prefix << " failed";
  }
  return test;
}
#endif


template<class Vector, class BlockingTag, class Coefficient, std::size_t dim, class MultiIndex>
Dune::TestSuite checkHierarchicVector(std::string shortName="")
{
  Dune::TestSuite test(shortName);

  using namespace Dune::TypeTree::Indices;
  using SizeInfo = HybridSizeInfoDummy<dim>;
  using SizePrefix = typename SizeInfo::SizePrefix;

  SizeInfo sizeInfo;

  // Create raw vector
  Vector x_raw;

  // Create wrapped vector
  auto x = Dune::Functions::hierarchicVectorWrapper(x_raw, BlockingTag{});

  // Resize wrapped vector using sizeInfo
  x.resize(sizeInfo);

  // Derive size information from vector
  test.require(x_raw.size() == sizeInfo(SizePrefix{}), "resize check")
    << "x_raw.size() is " << x_raw.size() << " but should be " << sizeInfo(SizePrefix{});

  test.require(x_raw[_0].size() == sizeInfo(SizePrefix{0}), "resize check")
    << "x_raw[_0].size() is " << x_raw[_0].size() << " but should be " << sizeInfo(SizePrefix{0});

  for (std::size_t i=0; i<sizeInfo({0}); ++i)
    test.require(x_raw[_0][i].size() == sizeInfo(SizePrefix{0,i}), "resize check")
      << "x_raw[_0][" << i << "].size() is " << x_raw[_0][i].size() << " but should be " << sizeInfo(SizePrefix{0,i});

  test.require(x_raw[_1].size() == sizeInfo(SizePrefix{1}), "resize check")
    << "x_raw[_1].size() is " << x_raw[_0].size() << " but should be " << sizeInfo(SizePrefix{1});


  // Assign values to each vector entry
  for (std::size_t i=0; i<x_raw[_0].size(); ++i)
    for (std::size_t j=0; j<x_raw[_0][i].size(); ++j)
      x[MultiIndex{{0, i, j}}] = 0+i+j;
  for (std::size_t i=0; i<x_raw[_1].size(); ++i)
    x[MultiIndex{{1, i}}] = 1+i;


  // Access vector entries via const reference
  const auto& x_const = x;
  for (std::size_t i=0; i<x_raw[_0].size(); ++i)
    for (std::size_t j=0; j<x_raw[_0][i].size(); ++j)
    {
      test.check(x_const[MultiIndex{{0, i, j}}] == Coefficient(0+i+j))
        << "x[{0," << i << "," << j << "}] contains wrong value";
    }
  for (std::size_t i=0; i<x_raw[_1].size(); ++i)
  {
    test.check(x_const[MultiIndex{{1, i}}] == Coefficient(1+i))
      << "x[{1," << i << "}] contains wrong value";
  }

  return test;
}


template<class Matrix, class RowBlocking, class ColBlocking, class Coefficient, std::size_t dim, class MultiIndex>
Dune::TestSuite checkHierarchicMatrix(std::string shortName="")
{
  Dune::TestSuite test(shortName);

  using namespace Dune::Functions;
  using namespace Dune::TypeTree::Indices;
  using SizeInfo = HybridSizeInfoDummy<dim>;
  using SizePrefix = typename SizeInfo::SizePrefix;

  SizeInfo sizeInfo;

  // Create raw vector
  Matrix mat_raw;

  // Create wrapped matrix
  auto mat = Dune::Functions::hierarchicMatrixWrapper(mat_raw, RowBlocking{}, ColBlocking{});

  // Resize wrapped matrix using sizeInfo
  mat.resize(sizeInfo, sizeInfo);

  std::cout << "size(mat) = " << hybridNumRows(mat_raw) << ", " << hybridNumCols(mat_raw) << "\n";
  std::cout << "size(mat[0][0]) = " << hybridNumRows(mat_raw[_0][_0]) << ", " << hybridNumCols(mat_raw[_0][_0]) << "\n";
  std::cout << "size(mat[0][1]) = " << hybridNumRows(mat_raw[_0][_1]) << ", " << hybridNumCols(mat_raw[_0][_1]) << "\n";
  std::cout << "size(mat[1][0]) = " << hybridNumRows(mat_raw[_1][_0]) << ", " << hybridNumCols(mat_raw[_1][_0]) << "\n";
  std::cout << "size(mat[1][1]) = " << hybridNumRows(mat_raw[_1][_1]) << ", " << hybridNumCols(mat_raw[_1][_1]) << "\n";

  // Derive size information from vector
  test.require(hybridNumRows(mat_raw) == sizeInfo(SizePrefix{}), "resize check")
    << "num_rows(mat_raw) is " << hybridNumRows(mat_raw) << " but should be " << sizeInfo(SizePrefix{});
  test.require(hybridNumCols(mat_raw) == sizeInfo(SizePrefix{}), "resize check")
    << "num_cols(mat_raw) is " << hybridNumCols(mat_raw) << " but should be " << sizeInfo(SizePrefix{});

  test.require(hybridNumRows(mat_raw[_0][_1]) == sizeInfo(SizePrefix{0}), "resize check")
    << "num_rows(mat_raw[_0][_1]) is " << hybridNumRows(mat_raw[_0][_1]) << " but should be " << sizeInfo(SizePrefix{0});

  test.require(hybridNumCols(mat_raw[_0][_1]) == sizeInfo(SizePrefix{1}), "resize check")
    << "num_cols(mat_raw[_0][_1]) is " << hybridNumCols(mat_raw[_0][_1]) << " but should be " << sizeInfo(SizePrefix{0});
/*
  for (std::size_t i=0; i<sizeInfo({0}); ++i)
    test.require(x_raw[_0][i].size() == sizeInfo(SizePrefix{0,i}), "resize check")
      << "x_raw[_0][" << i << "].size() is " << x_raw[_0][i].size() << " but should be " << sizeInfo(SizePrefix{0,i});

  test.require(x_raw[_1].size() == sizeInfo(SizePrefix{1}), "resize check")
    << "x_raw[_1].size() is " << x_raw[_0].size() << " but should be " << sizeInfo(SizePrefix{1});


  // Assign values to each vector entry
  for (std::size_t i=0; i<x_raw[_0].size(); ++i)
    for (std::size_t j=0; j<x_raw[_0][i].size(); ++j)
      x[MultiIndex{{0, i, j}}] = 0+i+j;
  for (std::size_t i=0; i<x_raw[_1].size(); ++i)
    x[MultiIndex{{1, i}}] = 1+i;


  // Access vector entries via const reference
  const auto& x_const = x;
  for (std::size_t i=0; i<x_raw[_0].size(); ++i)
    for (std::size_t j=0; j<x_raw[_0][i].size(); ++j)
    {
      test.check(x_const[MultiIndex{{0, i, j}}] == Coefficient(0+i+j))
        << "x[{0," << i << "," << j << "}] contains wrong value";
    }
  for (std::size_t i=0; i<x_raw[_1].size(); ++i)
  {
    test.check(x_const[MultiIndex{{1, i}}] == Coefficient(1+i))
      << "x[{1," << i << "}] contains wrong value";
  }*/

  return test;
}




int main (int argc, char *argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  using namespace Dune::Functions::Blocking;

  Dune::TestSuite test;
  {
    using VelocityVector = std::vector<std::vector<double>>;
    using PressureVector = std::vector<double>;
    using Vector = Dune::TupleVector<VelocityVector, PressureVector>;
    using Blocking = Blocked<Blocked<Flat>,Flat>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Blocking, double, 2, MultiIndex>("TV<V<V<double>>, V<double>>"));
  }

  {
    using VelocityVector = std::vector<Dune::BlockVector<Dune::FieldVector<double,1>>>;
    using PressureVector = std::vector<Dune::FieldVector<double,1>>;
    using Blocking = Blocked<Blocked<Flat>,Flat>;
    using Vector = Dune::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Blocking, double, 2, MultiIndex>("TV<V<BV<FV<double,1>>>, V<FV<doule,1>>>"));
  }

  {
    using VelocityVector = std::vector<std::vector<Dune::FieldVector<double,3>>>;
    using PressureVector = std::vector<Dune::FieldVector<double,3>>;
    using Blocking = Blocked<Blocked<Flat>,Flat>;
    using Vector = Dune::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Blocking, Dune::FieldVector<double,3>, 2, MultiIndex>("TV<V<V<FV<double,3>>>, V<FV<double,3>>>"));
  }

  {
    static const std::size_t dim = 5;
    using VelocityVector = std::vector<std::array<Dune::FieldVector<double,1>,dim>>;
    using PressureVector = std::vector<double>;
//     using Blocking = Blocked<Blocked<Blocked<Flat>>,Flat>; // or...
    using Blocking = Blocked<Blocked<LeafBlocked<1>>,Flat>;
    using Vector = Dune::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Blocking, double, dim, MultiIndex>("TV<V<A<FV<double,1>,5>>, V<double>>"));
  }

  {
    static const std::size_t dim = 5;
    using VelocityVector = Dune::BlockVector<Dune::FieldVector<double,dim>>;
    using PressureVector = Dune::BlockVector<Dune::FieldVector<double,1>>;
    using Blocking = Blocked<LeafBlocked<dim>,LeafBlocked<1>>;
    using Vector = Dune::MultiTypeBlockVector<VelocityVector, PressureVector>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Blocking, double, dim, MultiIndex>("MTBV<BV<FV<double,5>>, BV<FV<double,1>>>"));
  }

  { // does not yet work
#if 0
    static const std::size_t dim = 3;
    using VelocityVector = std::vector<Dune::MultiTypeBlockVector<Dune::FieldVector<double,1>, double, Dune::FieldVector<double,1>>>;
    using PressureVector = Dune::BlockVector<Dune::FieldVector<double,1>>;
    using Blocking = Blocked<LeafBlocked<3>,LeafBlocked<1>>;
    using Vector = Dune::MultiTypeBlockVector<VelocityVector, PressureVector>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Blocking, double, dim, MultiIndex>("MTBV<V<MTBV<FV<double,1>, double, FV<double,1>>>, BV<FV<double,1>>"));
#endif
  }

  // test with std::vector<bool>, i.e. no coefficient type bool
  {
    using VelocityVector = std::vector<std::vector<bool>>;
    using PressureVector = std::vector<bool>;
    using Vector = Dune::TupleVector<VelocityVector, PressureVector>;
    using Blocking = Blocked<Blocked<Flat>,Flat>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Blocking, bool, 2, MultiIndex>("TV<V<V<bool>>, V<bool>>"));
  }

  {
    using VelocityMatrix = Dune::Matrix<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>;
    using PressureMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
    using Matrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<VelocityMatrix, VelocityMatrix>, Dune::MultiTypeBlockVector<VelocityMatrix, PressureMatrix>>;
    using Blocking = Blocked<Blocked<Flat>,Flat>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicMatrix<Matrix, Blocking, Blocking, double, 2, MultiIndex>("TV<V<V<double>>, V<double>>"));
  }

  return test.exit();
}
// Error handling
catch (Dune::Exception e) {
  std::cout << e << std::endl;
  return 1;
}
