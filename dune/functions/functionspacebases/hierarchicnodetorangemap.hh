// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICNODETORANGEMAP_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICNODETORANGEMAP_HH


#include <utility>
#include <type_traits>

#include <dune/common/concept.hh>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/common/indexaccess.hh>

namespace Dune {
namespace Functions {



/**
 * \brief A simple node to range map using the nested tree indices
 *
 * This map directly usses the tree path entries of the given
 * node to access the nested container.
 *
 * If the container does not provide any operator[] access,
 * it is simply forwarded for all nodes.
 */
struct HierarchicNodeToRangeMap
{
  template<class Node, class TreePath, class Range,
    std::enable_if_t< models<Concept::HasIndexAccess, Range, Dune::index_constant<0>>(), int> = 0>
  decltype(auto) operator()(const Node& node, const TreePath& treePath, Range&& y) const
  {
    return resolveStaticMultiIndex(y, treePath);
  }

  template<class Node, class TreePath, class Range,
    std::enable_if_t<not models<Concept::HasIndexAccess, Range, Dune::index_constant<0>>(), int> = 0>
  decltype(auto) operator()(const Node& node, const TreePath& treePath, Range&& y) const
  {
    return std::forward<Range>(y);
  }
};



} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICNODETORANGEMAP_HH
