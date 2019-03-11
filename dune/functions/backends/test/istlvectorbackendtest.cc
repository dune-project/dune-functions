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

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dune/typetree/utility.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/backends/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>


using namespace Dune;



/**
 * \brief A Dummy global basis
 *
 * This is a mock class providing non-uniform size information.
 * It's non-uniform in the sense that not all multi-indices
 * have the same length.
 */
template<std::size_t dim>
class GlobalBasisMoc
{
public:
    using size_type = std::size_t;
    using SizePrefix = Dune::ReservedVector<std::size_t, 3>;
    using MultiIndex = Dune::ReservedVector<std::size_t, 3>;

    /**
     * \brief Construct from basis
     */
    GlobalBasisMoc()
    {}

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
      return 0;
    }

    operator size_type () const
    {
        return 23*dim+42;
    }

};


template<class Vector, class Coefficient, std::size_t dim, class MultiIndex>
Dune::TestSuite checkISTLVectorBackend(std::string shortName="")
{
  Dune::TestSuite test(shortName);

  using namespace Dune::TypeTree::Indices;
  using Basis = GlobalBasisMoc<dim>;
  using SizePrefix = typename Basis::SizePrefix;

  Basis basis;

  // Create raw vector
  Vector x_raw;

  // Create wrapped vector
  auto x = Dune::Functions::istlVectorBackend(x_raw);

  test.require(Dune::models<Dune::Functions::Concept::VectorBackend<Basis>, decltype(x)>(), "VectorBackend concept check")
    << "Object returned by istlVectorBackend() does not model the VectorBackend concept";

  // Resize wrapped vector using basis
  x.resize(basis);

  // Derive size information from vector
  test.require(x_raw.size() == basis.size(SizePrefix{}), "resize check")
    << "x_raw.size() is " << x_raw.size() << " but should be " << basis.size(SizePrefix{});

  test.require(x_raw[_0].size() == basis.size(SizePrefix{0}), "resize check")
    << "x_raw[_0].size() is " << x_raw[_0].size() << " but should be " << basis.size(SizePrefix{0});

  for (std::size_t i=0; i<basis.size({0}); ++i)
    test.require(x_raw[_0][i].size() == basis.size(SizePrefix{0,i}), "resize check")
      << "x_raw[_0][" << i << "].size() is " << x_raw[_0][i].size() << " but should be " << basis.size(SizePrefix{0,i});

  test.require(x_raw[_1].size() == basis.size(SizePrefix{1}), "resize check")
    << "x_raw[_1].size() is " << x_raw[_0].size() << " but should be " << basis.size(SizePrefix{1});


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
    using Vector = Dune::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkISTLVectorBackend<Vector, Coefficient, 2, MultiIndex>("TV<V<V<double>>, V<double>>"));
  }

  {
    using VelocityVector = std::vector<Dune::BlockVector<Dune::FieldVector<double,1>>>;
    using PressureVector = std::vector<Dune::FieldVector<double,1>>;
    using Coefficient = double;
    using Vector = Dune::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkISTLVectorBackend<Vector, Coefficient, 2, MultiIndex>("TV<V<BV<FV<double,1>>>, V<FV<doule,1>>>"));
  }

  {
    static const std::size_t dim = 5;
    using VelocityVector = std::vector<std::array<Dune::FieldVector<double,1>,dim>>;
    using PressureVector = std::vector<double>;
    using Coefficient = double;
    using Vector = Dune::TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkISTLVectorBackend<Vector, Coefficient, dim, MultiIndex>("TV<V<A<FV<double,1>,5>>, V<double>>"));
  }

  {
    static const std::size_t dim = 5;
    using VelocityVector = Dune::BlockVector<Dune::FieldVector<double,dim>>;
    using PressureVector = Dune::BlockVector<Dune::FieldVector<double,1>>;
    using Coefficient = double;
    using Vector = Dune::MultiTypeBlockVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkISTLVectorBackend<Vector, Coefficient, dim, MultiIndex>("MTBV<BV<FV<double,5>>, BV<FV<double,1>>>"));
  }

  {
    static const std::size_t dim = 3;
    using VelocityVector = std::vector<Dune::MultiTypeBlockVector<Dune::FieldVector<double,1>, double, Dune::FieldVector<double,1>>>;
    using PressureVector = Dune::BlockVector<Dune::FieldVector<double,1>>;
    using Coefficient = double;
    using Vector = Dune::MultiTypeBlockVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    test.subTest(checkISTLVectorBackend<Vector, Coefficient, dim, MultiIndex>("MTBV<V<MTBV<FV<double,1>, double, FV<double,1>>>, BV<FV<double,1>>"));
  }


  return test.exit();
}
// Error handling
catch (Exception& e) {
  std::cout << e.what() << std::endl;
  return 1;
}
