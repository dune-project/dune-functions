// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>
#include <array>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/classname.hh>
#include <dune/common/fvector.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dune/typetree/utility.hh>

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




/**
 * \brief A simple multi-type container
 */
template<class... T>
class TupleVector : public std::tuple<T...>
{
  using Base = std::tuple<T...>;

public:

  template<class... TT>
  constexpr TupleVector(TT&&... tt) :
    Base(std::forward<TT>(tt)...)
  {}

  constexpr TupleVector()
  {}

  template<std::size_t i>
  auto operator[](const Dune::TypeTree::index_constant<i>&) const
    ->decltype(std::get<i>(*this))
  {
    return std::get<i>(*this);
  }

  template<std::size_t i>
  auto operator[](const Dune::TypeTree::index_constant<i>&)
    ->decltype(std::get<i>(*this))
  {
    return std::get<i>(*this);
  }

  static constexpr std::size_t size()
  {
    return std::tuple_size<Base>::value;
  }

};


class TestResult
{

  class TestStream : public std::ostringstream
  {
  public:

    TestStream(bool condition, bool required) :
      condition_(condition),
      required_(required)
    {}

    TestStream(const TestStream& other) :
      condition_(other.condition_),
      required_(other.required_)
    {
      (*this) << other.str();
      other.condition_ = true;
    }

    ~TestStream()
    {
      if (not condition_)
      {
        if (required_)
        {
          std::cout << "Required check failed : " << this->str() << std::endl;
          DUNE_THROW(Dune::Exception, "Required check failed : " << this->str());
        }
        else
          std::cout << "Check failed : " << this->str() << std::endl;
      }
    }

  private:
    mutable bool condition_;
    bool required_;
  };

public:

  TestResult() :
    result_(true)
  {}

  TestStream check(bool condition)
  {
    result_ &= condition;
    return TestStream(condition, false);
  }

  TestStream require(bool condition)
  {
    result_ &= condition;
    return TestStream(condition, true);
  }

  operator const bool& () const
  {
    return result_;
  }

  operator bool& ()
  {
    return result_;
  }

private:
  bool result_;
};



template<class Vector, class Coefficient, std::size_t dim, class MultiIndex>
bool checkHierarchicVector()
{
  TestResult result;

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
  result.require(x_raw.size() == sizeInfo(SizePrefix{}))
    << "x_raw.size() is " << x_raw.size() << " but should be " << sizeInfo(SizePrefix{});

  result.require(x_raw[_0].size() == sizeInfo(SizePrefix{0}))
    << "x_raw[_0].size() is " << x_raw[_0].size() << " but should be " << sizeInfo(SizePrefix{0});

  for (std::size_t i=0; i<sizeInfo({0}); ++i)
    result.require(x_raw[_0][i].size() == sizeInfo(SizePrefix{0,i}))
      << "x_raw[_0][" << i << "].size() is " << x_raw[_0][i].size() << " but should be " << sizeInfo(SizePrefix{0,i});

  result.require(x_raw[_1].size() == sizeInfo(SizePrefix{1}))
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
      std::cout << "x[{0," << i << "," << j << "}] = " << x_const[MultiIndex{{0, i, j}}] << std::endl;
//      result.check(x[MultiIndex{{0, i, j}}] == 0.0+i+j)
//        << "x[{0," << i << "," << j << "}] contains wrong value";
    }
  for (std::size_t i=0; i<x_raw[_1].size(); ++i)
  {
    std::cout << "x[{1," << i << "}] = " << x_const[MultiIndex{{1, i}}] << std::endl;
//    result.check(x[MultiIndex{{1, i}}] == 1.0+i)
//      << "x[{1," << i << "}] contains wrong value";
  }
  return result;
}


int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////


  TestResult result;
  {
    using VelocityVector = std::vector<std::vector<double>>;
    using PressureVector = std::vector<double>;
    using Coefficient = double;
    using Vector = TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    result.check(checkHierarchicVector<Vector, Coefficient, 2, MultiIndex>())
      << "Test with " << Dune::className<Vector>() << " failed";
  }

  {
    using VelocityVector = std::vector<Dune::BlockVector<Dune::FieldVector<double,1>>>;
    using PressureVector = std::vector<Dune::FieldVector<double,1>>;
    using Coefficient = double;
    using Vector = TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    result.check(checkHierarchicVector<Vector, Coefficient, 2, MultiIndex>())
      << "Test with " << Dune::className<Vector>() << " failed";
  }

  {
    using VelocityVector = std::vector<std::vector<Dune::FieldVector<double,3>>>;
    using PressureVector = std::vector<Dune::FieldVector<double,3>>;
    using Coefficient = Dune::FieldVector<double,3>;
    using Vector = TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    result.check(checkHierarchicVector<Vector, Coefficient, 2, MultiIndex>())
      << "Test with " << Dune::className<Vector>() << " failed";
  }

  {
    static const std::size_t dim = 5;
    using VelocityVector = std::vector<std::array<Dune::FieldVector<double,1>,dim>>;
    using PressureVector = std::vector<double>;
    using Coefficient = double;
    using Vector = TupleVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    result.check(checkHierarchicVector<Vector, Coefficient, dim, MultiIndex>())
      << "Test with " << Dune::className<Vector>() << " failed";
  }

  {
    static const std::size_t dim = 5;
    using VelocityVector = Dune::BlockVector<Dune::FieldVector<double,dim>>;
    using PressureVector = Dune::BlockVector<Dune::FieldVector<double,1>>;
    using Coefficient = double;
    using Vector = Dune::MultiTypeBlockVector<VelocityVector, PressureVector>;
    using MultiIndex = ReservedVector<std::size_t, 3>;
    result.check(checkHierarchicVector<Vector, Coefficient, dim, MultiIndex>())
      << "Test with " << Dune::className<Vector>() << " failed";
  }

  if (not result)
    std::cout << "Test failed" << std::endl;

  return result ? 0 : 1;

}
// Error handling
catch (Exception e) {
  std::cout << e << std::endl;
}
