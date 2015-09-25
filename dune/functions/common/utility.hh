// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_UTILITY_HH
#define DUNE_FUNCTIONS_COMMON_UTILITY_HH


#include <utility>
#include <type_traits>

#include <dune/common/std/utility.hh>
#include <dune/typetree/utility.hh>

namespace Dune {
namespace Functions {



template<class F, class size_type, size_type firstValue, size_type secondValue, size_type... otherValues, class... Args>
auto forwardAsStaticInteger(Dune::Std::integer_sequence<size_type, firstValue, secondValue, otherValues...> values, const size_type i, F&& f, Args&&... args)
  ->decltype(f(std::integral_constant<size_type, firstValue>(), std::forward<Args>(args)...))
{
  if (i==firstValue)
    return f(std::integral_constant<size_type, firstValue>(), std::forward<Args>(args)...);
  else
    return forwardAsStaticInteger(Dune::Std::integer_sequence<size_type, secondValue, otherValues...>(), i, std::forward<F>(f), std::forward<Args>(args)...);
}

template<class F, class size_type, size_type firstValue, class... Args>
auto forwardAsStaticInteger(Dune::Std::integer_sequence<size_type, firstValue> values, const size_type& i, F&& f, Args&&... args)
  ->decltype(f(std::integral_constant<size_type, firstValue>(), std::forward<Args>(args)...))
{
  return f(std::integral_constant<size_type, firstValue>(), std::forward<Args>(args)...);
}

/**
 * \brief Transform dynamic index to static index_constant
 *
 * This will call the given function with index_constant<i>
 * where i is the dynamically provided index.
 *
 * To achieve this the conditon i==ii is checked subsequently
 * for al static indices ii in the range 0,...,(end-1). In order
 * to be able to compile this we require for all ii in this range
 * that f(index_constant<ii>()) is well-formed and that the result
 * type of it can be converted to the result type of f(index_constant<0>()).
 * If i is not in this range, the returned value is f(index_constant<n-1>())
 *
 * \param i Dynamic index
 * \param f Function to call (e.g., a generic lambda)
 * \param args Additional arguments for f
 *
 * \returns f(index_constant<i>(), args...)
 */
template<std::size_t end, class F, class size_type, class... Args>
auto forwardAsStaticIndex(const size_type& i, F&& f, Args&&... args)
  ->decltype(f(Dune::TypeTree::Indices::_0, std::forward<Args>(args)...))
{
  return forwardAsStaticInteger(Dune::TypeTree::Std::make_index_sequence<end>(), i, std::forward<F>(f), std::forward<Args>(args)...);
}


} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_COMMON_UTILITY_HH
