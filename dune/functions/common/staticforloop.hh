// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_STATICFORLOOP_HH
#define DUNE_FUNCTIONS_COMMON_STATICFORLOOP_HH


#include <dune/common/concept.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/concepts.hh>


namespace Dune {
namespace Functions {

namespace Imp {

template<class ST, ST begin, ST end>
struct StaticFindInRange
{
  template<class F, class...Args>
  static void apply(F&& f, Args&&... args)
  {
    if (f(std::integral_constant<ST, begin>(), std::forward<Args>(args)...))
      return;
    StaticFindInRange<ST, begin+1, end>::apply(std::forward<F>(f), std::forward<Args>(args)...);
  }
};

template<class ST, ST end>
struct StaticFindInRange<ST, end, end>
{
  template<class F, class...Args>
  static void apply(F&& f, Args&&...)
  {}
};

} //end namespace Imp



/**
 * \brief Static find loop
 *
 * \ingroup Utility
 *
 * Run static for-loop from 'begin' to 'end-1' with functor.
 * The functor is called with \p TypeTree::index_constant<i>
 * as first argument. All other arguments of this method
 * are forwarded to the functor. If the functor returns
 * true the loop is terminated.
 *
 * \todo Should this be just the StaticForLoop?
 */
template<std::size_t begin_t, std::size_t end_t, class F, class... Args>
void staticFindInRange(F&& f, Args&&... args)
{
  Imp::StaticFindInRange<std::size_t, begin_t, end_t>::apply(std::forward<F>(f), std::forward<Args>(args)...);
}


} // namespace Dune::Functions
} // namespace Dune



#endif //DUNE_FUNCTIONS_COMMON_STATICFORLOOP_HH
