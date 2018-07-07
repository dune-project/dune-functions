// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GEOMETRYTYPEPROVIDER_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GEOMETRYTYPEPROVIDER_HH

#include <dune/common/std/type_traits.hh>
#include <dune/common/indices.hh>
#include <dune/common/keywords.hh>
#include <dune/geometry/type.hh>

#include <dune/functions/common/staticgeometrytype.hh>



namespace Dune {
namespace Functions {



struct MixedGeometryTypeProvider
{
  template<class GridView, class Entity>
  static constexpr auto type(const Entity& entity)
  {
    return entity.type();
  }

  template<class GridView>
  using Type = GeometryType;
};

struct SimplexGeometryTypeProvider
{
  template<class GridView, class Entity>
  static constexpr auto type(const Entity& entity)
  {
    return StaticGeometryTypes::Simplex<GridView::dimension>();
  }

  template<class GridView>
  using Type = StaticGeometryTypes::Simplex<GridView::dimension>;
};

struct CubeGeometryTypeProvider
{
  template<class GridView, class Entity>
  static constexpr auto type(const Entity& entity)
  {
    return StaticGeometryTypes::Cube<GridView::dimension>();
  }

  template<class GridView>
  using Type = StaticGeometryTypes::Cube<GridView::dimension>;
};

struct AutoGeometryTypeProvider
{
  template<class GridView, class Entity,
    std::enable_if_t<Dune::Capabilities::hasSingleGeometryType<typename GridView::Grid>::v, int> = 0>
  static constexpr auto type(const Entity& entity)
  {
    constexpr unsigned int topologyId = Dune::Capabilities::hasSingleGeometryType<typename GridView::Grid>::topologyId;
    constexpr unsigned int normalizedTopologyId = (topologyId >> 1) << 1;
    return StaticGeometryType<normalizedTopologyId, GridView::dimension, false>{};
  }

  template<class GridView, class Entity,
    std::enable_if_t<not Dune::Capabilities::hasSingleGeometryType<typename GridView::Grid>::v, int> = 0>
  static constexpr auto type(const Entity& entity)
  {
    return entity.type();
  }

  template<class GridView>
  using Type = decltype(type<GridView>(std::declval<typename GridView::template Codim<0>::Entity>()));
};



} // namespace Functions
} // namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GEOMETRYTYPEPROVIDER_HH
