// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_STATICFORLOOP_HH
#define DUNE_FUNCTIONS_COMMON_STATICFORLOOP_HH


#include <dune/functions/common/concept.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/concepts.hh>


namespace Dune {
namespace Functions {



template<class ST, ST begin, ST end>
struct StaticForLoop
{
  template<class F, class...Args>
  static void apply(F&& f, Args&&... args)
  {
    f(std::integral_constant<ST, begin>(), std::forward<Args>(args)...);
    StaticForLoop<ST, begin+1, end>::apply(std::forward<F>(f), std::forward<Args>(args)...);
  }
};

template<class ST, ST end>
struct StaticForLoop<ST, end, end>
{
  template<class F, class...Args>
  static void apply(F&& f, Args&&...)
  {}
};

/**
 * \brief Static for loop
 *
 * Run static for-loop from 'begin' to 'end-1' with functor.
 * The functor is called with TypeTree::index_constant<i>
 * as first argument. All other arguments of this method
 * are forwarded to the functor.
 */
template<std::size_t begin_t, std::size_t end_t, class F, class... Args>
void staticForLoop(F&& f, Args&&... args)
{
  StaticForLoop<std::size_t, begin_t, end_t>::apply(std::forward<F>(f), std::forward<Args>(args)...);
}



} // namespace Dune::Functions
} // namespace Dune



#endif //DUNE_FUNCTIONS_COMMON_STATICFORLOOP_HH
