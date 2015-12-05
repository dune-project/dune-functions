// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>
#include <array>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/classname.hh>
#include <dune/common/fvector.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dune/typetree/utility.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/tuplevector.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>


using namespace Dune;



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


template<class Vector, class Coefficient, std::size_t dim, class MultiIndex>
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
  Dune::Functions::HierarchicVectorWrapper<Vector, Coefficient> x(x_raw);

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


int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////


  Dune::TestSuite test;
  {
    using VelocityVector = std::vector<std::vector<double>>;
    using PressureVector = std::vector<double>;
    using Coefficient = double;
    using Vector = Functions::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Coefficient, 2, MultiIndex>("TV<V<V<double>>, V<double>>"));
  }

  {
    using VelocityVector = std::vector<Dune::BlockVector<Dune::FieldVector<double,1>>>;
    using PressureVector = std::vector<Dune::FieldVector<double,1>>;
    using Coefficient = double;
    using Vector = Functions::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Coefficient, 2, MultiIndex>("TV<V<BV<FV<double,1>>>, V<FV<doule,1>>>"));
  }

  {
    using VelocityVector = std::vector<std::vector<Dune::FieldVector<double,3>>>;
    using PressureVector = std::vector<Dune::FieldVector<double,3>>;
    using Coefficient = Dune::FieldVector<double,3>;
    using Vector = Functions::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Coefficient, 2, MultiIndex>("TV<V<V<FV<double,3>>>, V<FV<double,3>>>"));
  }

  {
    static const std::size_t dim = 5;
    using VelocityVector = std::vector<std::array<Dune::FieldVector<double,1>,dim>>;
    using PressureVector = std::vector<double>;
    using Coefficient = double;
    using Vector = Functions::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Coefficient, dim, MultiIndex>("TV<V<A<FV<double,1>,5>>, V<double>>"));
  }

  {
    static const std::size_t dim = 5;
    using VelocityVector = Dune::BlockVector<Dune::FieldVector<double,dim>>;
    using PressureVector = Dune::BlockVector<Dune::FieldVector<double,1>>;
    using Coefficient = double;
    using Vector = Dune::MultiTypeBlockVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Coefficient, dim, MultiIndex>("MTBV<BV<FV<double,5>>, BV<FV<double,1>>>"));
  }

  {
    static const std::size_t dim = 3;
    using VelocityVector = std::vector<Dune::MultiTypeBlockVector<Dune::FieldVector<double,1>, double, Dune::FieldVector<double,1>>>;
    using PressureVector = Dune::BlockVector<Dune::FieldVector<double,1>>;
    using Coefficient = double;
    using Vector = Dune::MultiTypeBlockVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkHierarchicVector<Vector, Coefficient, dim, MultiIndex>("MTBV<V<MTBV<FV<double,1>, double, FV<double,1>>>, BV<FV<double,1>>"));
  }


  return test.exit();
}
// Error handling
catch (Exception e) {
  std::cout << e << std::endl;
  return 1;
}
