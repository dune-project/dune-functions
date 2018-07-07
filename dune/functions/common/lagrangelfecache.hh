// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_LAGRANGELFECACHE_HH
#define DUNE_FUNCTIONS_COMMON_LAGRANGELFECACHE_HH

#include <map>

#include <dune/common/tuplevector.hh>
#include <dune/common/overloadset.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/prismp1.hh>
#include <dune/localfunctions/lagrange/prismp2.hh>
#include <dune/localfunctions/lagrange/pyramidp1.hh>
#include <dune/localfunctions/lagrange/pyramidp2.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/functions/common/staticgeometrytype.hh>



namespace Dune {
namespace Functions {



namespace Imp {

  // This class is used as fallback if there's no matching Lagrange
  // finite element implementation. By having the required static
  // geometry type and order as template parameter, the compiler
  // will print out this information in the error message in case
  // this dummy is used. This is helpful when debugging the reason.
  template<class SGT, std::size_t order>
  struct DummyLocalFiniteElement {};



  // Create Lagrange finite element for given geometry type, ctypes and order
  template<class D, class R, unsigned int order, unsigned int topologyId, std::size_t dim, bool none>
  constexpr auto lagrangeFiniteElement(const StaticGeometryType<topologyId, dim, none>& type)
  {
    // Using orderedOverload allows to avoid explicit SFINAE
    // conditions, because the first match is used which
    // avoids ambiguity by construction.
    return Dune::orderedOverload(
      [](Dune::index_constant<0> o, auto t) {
        return P0LocalFiniteElement<D,R,dim>(t);
      },
      [](auto o, const StaticGeometryTypes::Simplex<dim>& t) {
        return PkLocalFiniteElement<D,R,dim,o>();
      },
      [](auto o, const StaticGeometryTypes::Cube<dim>& t) {
        return QkLocalFiniteElement<D,R,dim,o>();
      },
      [](Dune::index_constant<1> o, StaticGeometryTypes::Prism t) {
        return PrismP1LocalFiniteElement<D,R>();
      },
      [](Dune::index_constant<2> o, StaticGeometryTypes::Prism t) {
        return PrismP2LocalFiniteElement<D,R>();
      },
      [](Dune::index_constant<1> o, StaticGeometryTypes::Pyramid t) {
        return PyramidP1LocalFiniteElement<D,R>();
      },
      [](Dune::index_constant<2> o, StaticGeometryTypes::Pyramid t) {
        return PyramidP2LocalFiniteElement<D,R>();
      },
      [](auto o, auto t) {
        return Imp::DummyLocalFiniteElement<decltype(t), o>();
      }
      )(Dune::index_constant<order>(), type);
  }


} // namespace Imp



/** \brief A cache that stores all available Pk/Qk like local finite elements for the given dimension and order
 *
 * In contrast to the PQkLocalFiniteElementCache in dune-localfunctions this
 * class also provides access using StaticGeometryType which will directly
 * provide the finite element implementation without wrapping it into the
 * polymorphic interface.
 *
 * An interface for dealing with different vertex orders is currently missing.
 *
 * \tparam D Type used for domain coordinates
 * \tparam R Type used for shape function values
 * \tparam dim Element dimension
 * \tparam order Element order
 */
template<class D, class R, int dim, int order>
class LagrangeFiniteElementCache
{

  static constexpr auto createStaticCache()
  {
    auto gtIndices = std::make_index_sequence<LocalHybridGeometryTypeIndex::size(dim)>();
    return unpackIntegerSequence([&](auto... index) {
      return Dune::makeTupleVector(
          Imp::lagrangeFiniteElement<D, R, order>(LocalHybridGeometryTypeIndex::type<dim, index>())...);
    }, gtIndices);
  }

  using StaticCache = decltype(createStaticCache());

  using DynamicCache = Dune::PQkLocalFiniteElementCache<D, R, dim, order>;

public:

  LagrangeFiniteElementCache() :
    dynamicCache_(),
    staticCache_(createStaticCache())
  {}

  LagrangeFiniteElementCache(const LagrangeFiniteElementCache& other) = default;

  template<unsigned int topologyId>
  const auto& get(const StaticGeometryType<topologyId, dim, false>& gt) const
  {
    return staticCache_[LocalHybridGeometryTypeIndex::index(gt)];
  }

  const auto& get(const GeometryType& gt) const
  {
    return dynamicCache_.get(gt);
  }

  template<class GT>
  using FiniteElement = std::decay_t<decltype(std::declval<LagrangeFiniteElementCache>().get(std::declval<GT>()))>;

private:

  DynamicCache dynamicCache_;
  StaticCache staticCache_;
};



} // namespace Functions
} // namespace Dune




#endif // DUNE_FUNCTIONS_COMMON_LAGRANGELFECACHE_HH
