// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

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



int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  static const int dim = 2;

  using namespace Dune::TypeTree::Indices;

//  using VelocityVector = std::vector<std::array<double,dim>>;
//  using PressureVector = std::vector<double>;

  using VelocityVector = std::vector<std::vector<double>>;
  using PressureVector = std::vector<double>;

  using Vector = TupleVector<VelocityVector, PressureVector>;

  Vector x_raw;

  auto x = Dune::Functions::HierarchicVectorWrapper<decltype(x_raw),double>(x_raw);
//  auto x = Dune::Functions::hierarchicVector(x_raw);

  HybridSizeInfoDummy<dim> sizeInfo;

  x.resize(sizeInfo);

  using MultiIndex = ReservedVector<std::size_t, 3>;
//  using MultiIndex = std::array<std::size_t, 3>;


  std::cout << "size()       " << x_raw.size() << std::endl;
  std::cout << "size({0})    " << x_raw[_0].size() << std::endl;
  for (std::size_t i=0; i<sizeInfo({0}); ++i)
    std::cout << "size({0," << i << "})  " << x_raw[_0][i].size() << std::endl;
  std::cout << "size({1})    " << x_raw[_1].size() << std::endl;


  for (std::size_t i=0; i<x_raw[_0].size(); ++i)
    for (std::size_t j=0; j<x_raw[_0][i].size(); ++j)
      x[MultiIndex{{0, i, j}}] = 0+i+j;

  for (std::size_t i=0; i<x_raw[_1].size(); ++i)
    x[MultiIndex{{1, i}}] = 1+i;



  const auto& x_const = x;

  for (std::size_t i=0; i<x_raw[_0].size(); ++i)
    for (std::size_t j=0; j<x_raw[_0][i].size(); ++j)
      std::cout << "x[{0," << i << "," << j << "}] = " << x_const[MultiIndex{{0, i, j}}] << std::endl;

  for (std::size_t i=0; i<x_raw[_1].size(); ++i)
      std::cout << "x[{1," << i << "}] = " << x_const[MultiIndex{{1, i}}] << std::endl;

}
// Error handling
catch (Exception e) {
  std::cout << e << std::endl;
}
