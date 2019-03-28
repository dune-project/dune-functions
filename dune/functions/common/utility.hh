// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_UTILITY_HH
#define DUNE_FUNCTIONS_COMMON_UTILITY_HH


#include <utility>
#include <type_traits>

#include <dune/common/overloadset.hh>
#include <dune/typetree/utility.hh>

#include <dune/functions/common/functionconcepts.hh>

namespace Dune {
namespace Functions {



template<class F, class size_type, size_type firstValue, class... Args>
auto forwardAsStaticInteger(std::integer_sequence<size_type, firstValue> values, const size_type& i, F&& f, Args&&... args)
  ->decltype(f(std::integral_constant<size_type, firstValue>(), std::forward<Args>(args)...))
{
  return f(std::integral_constant<size_type, firstValue>(), std::forward<Args>(args)...);
}

template<class F, class size_type, size_type firstValue, size_type secondValue, size_type... otherValues, class... Args>
auto forwardAsStaticInteger(std::integer_sequence<size_type, firstValue, secondValue, otherValues...> values, const size_type i, F&& f, Args&&... args)
  ->decltype(f(std::integral_constant<size_type, firstValue>(), std::forward<Args>(args)...))
{
  if (i==firstValue)
    return f(std::integral_constant<size_type, firstValue>(), std::forward<Args>(args)...);
  return forwardAsStaticInteger(std::integer_sequence<size_type, secondValue, otherValues...>(), i, std::forward<F>(f), std::forward<Args>(args)...);
}



/**
 * \brief Transform dynamic index to static index_constant
 *
 * \ingroup Utility
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
  return forwardAsStaticInteger(std::make_index_sequence<end>{}, i, std::forward<F>(f), std::forward<Args>(args)...);
}



namespace Imp {

  template<template<class...> class T, class List>
  struct ExpandTupleHelper
  {};

  template<template<class...> class T, template<class...> class ListType, class... Args>
  struct ExpandTupleHelper<T, ListType<Args...>>
  {
    using Type = T<Args...>;
  };

} // end namespace Imp

/**
 * \brief Expand tuple arguments as template arguments
 *
 * \ingroup Utility
 *
 * This template alias refers to T<Args...> if
 * ArgTuple is a std::tuple<Args...>.
 *
 * \tparam T A variadic template
 * \tparam ArgTuple A tuple of types
 */
template<template<class...> class T, class ArgTuple>
using ExpandTuple = typename Imp::ExpandTupleHelper<T, ArgTuple>::Type;



namespace Imp {

  template<template<class...> class T, class... Tuple>
  struct TransformTupleHelper
  {};

  template<template<class...> class T, class... Args1>
  struct TransformTupleHelper<T, typename std::tuple<Args1...>>
  {
    using Type = std::tuple<T<Args1>...>;
  };

  template<template<class...> class T, class... Args1, class... Args2>
  struct TransformTupleHelper<T, typename std::tuple<Args1...>, typename std::tuple<Args2...>>
  {
    using Type = std::tuple<T<Args1, Args2>...>;
  };

} // end namespace Imp

/**
 * \brief Transform tuple types argument using type-functor
 *
 * \ingroup Utility
 *
 * This is a template alias for a tuple whose i-th type
 * is given by F<T1i,...,TMi> where T1i,...,TMi are the
 * i-th types of the 1,...,M-th tuple of the given tuple
 * list Tuples. Currently only M=1,2 are supported.
 * \tparam F A template alias mapping 1,...,sizeof...(ArgTuple) types to a new one
 * \tparam Tuples A list of tuples
 */
template<template<class...> class F, class... Tuples>
using TransformTuple = typename Imp::TransformTupleHelper<F, Tuples...>::Type;



namespace Imp {

  template<class F, class... T, std::size_t... k>
  auto transformTupleHelper(F&& f, const std::tuple<T...>& tuple, std::index_sequence<k...>)
    -> decltype(std::make_tuple(f(std::get<k>(tuple))...))
  {
    return std::make_tuple(f(std::get<k>(tuple))...);
  }

  template<class F, class... T1, class...T2, std::size_t... k>
  auto transformTupleHelper(F&& f, const std::tuple<T1...>& tuple1, const std::tuple<T2...>& tuple2, std::index_sequence<k...>)
    -> decltype(std::make_tuple(f(std::get<k>(tuple1), std::get<k>(tuple2))...))
  {
    return std::make_tuple(f(std::get<k>(tuple1), std::get<k>(tuple2))...);
  }

} // end namespace Imp

/**
 * \brief Transform tuple value using a functor
 *
 * \ingroup Utility
 *
 * This will apply the given functor to all values in
 * given tuple and return the results in a new tuple.
 *
 * \param F A functor defined for all tuple entries
 * \param tuple The tuple to transform
 */
template<class F, class... T>
auto transformTuple(F&& f, const std::tuple<T...>& tuple)
  -> decltype(Imp::transformTupleHelper(std::forward<F>(f), tuple, std::index_sequence_for<T...>{}))
{
  return Imp::transformTupleHelper(std::forward<F>(f), tuple, std::index_sequence_for<T...>{});
}

/**
 * \brief Transform tuple value using a binary functor
 *
 * \ingroup Utility
 *
 * This will apply the given functor to the each corresponding
 * pair of values in the given tuples and return the results
 * in a new tuple.
 *
 * \param F A functor defined for all tuple entries
 * \param tuple1 The tuple containing values for the first parameter
 * \param tuple2 The tuple containing values for the second parameter
 */
template<class F, class... T1, class... T2>
auto transformTuple(F&& f, const std::tuple<T1...>& tuple1, const std::tuple<T2...>& tuple2)
  -> decltype(Imp::transformTupleHelper(std::forward<F>(f), tuple1, tuple2, std::index_sequence_for<T1...>{}))
{
  return Imp::transformTupleHelper(std::forward<F>(f), tuple1, tuple2, std::index_sequence_for<T1...>{});
}



namespace Imp {

  template<class IntegerSequence>
  struct IntegerSequenceTupleHelper
  {};

  template<class I, I... k>
  struct IntegerSequenceTupleHelper<std::integer_sequence<I, k...>>
  {
    using Type = std::tuple<std::integral_constant<I, k>...>;
  };

} // end namespace Imp

/**
 * \brief Transform integer_sequence<I,k...> to tuple<integral_constant<I,k>...>
 */
template<class IntegerSequence>
using IntegerSequenceTuple= typename Imp::IntegerSequenceTupleHelper<IntegerSequence>::Type;



/**
 * \brief Get last entry of type list
 *
 * \ingroup Utility
 */
template<class... T>
struct LastType
{
  using type = typename std::tuple_element<sizeof...(T)-1, std::tuple<T...>>::type;
};



namespace Imp {

template<class T, class I>
struct RotateHelper;

template<class... T, std::size_t... I>
struct RotateHelper<std::tuple<T...>, std::index_sequence<I...> >
{
  using type = typename std::tuple<typename LastType<T...>::type, typename std::tuple_element<I,std::tuple<T...>>::type...>;
};

} // end namespace Imp


/**
 * \brief Rotate type list by one, such that last entry is moved to first position
 *
 * \ingroup Utility
 *
 * The rotated type list is exported as tuple
 */
template<class... T>
struct RotateTuple
{
  using type = typename Imp::RotateHelper<std::tuple<T...>, std::make_index_sequence<sizeof...(T)-1>>::type;
};



/**
 * \brief Create a predicate for checking validity of expressions
 *
 * \param f A function involving the expression to check.
 *
 * This returns a function object that allows to check if the
 * expression encoded in f is valid for the given arguments.
 * To be precise it checks if f can be called using the given arguments.
 * This can be used inthe following way: To generate a check if the
 * expression x(a,b) is valid for some a and b use:
 *
 \code{.cpp}
 auto xIsValid = callableCheck([](auto&& a, auto&& b) -> void_t<decltype(x(a,b))> {});
 if (xIsValid(a,b))
   ...
 \endcode
 *
 * Notice that the given function f is stored by value.
 *
 * \ingroup Utility
 */
template<class Expression>
auto callableCheck(Expression f)
{
  return [f](auto&&... args){
    return Functions::Concept::isCallable(f, std::forward<decltype(args)>(args)...);
  };
}



/**
 * \brief Negate given predicate
 *
 * \param f A predicate function to negate
 *
 * This returns a function havin the same parameters as
 * f, but negating the result. Negation here means that
 * std::true_type is converted to std::false_type are
 * vice verse, while other return values are converted to
 * bool values and then the negated value is returned as bool, too.
 *
 * Notice that the given function f is stored by value.
 *
 * \ingroup Utility
 */
template<class Check>
auto negatePredicate(Check check)
{
  return [check](auto&&... args){
    auto negate = overload(
        [](std::true_type) { return std::false_type{};},
        [](std::false_type) { return std::true_type{};},
        [](bool v) { return not v;});
    return negate(check(std::forward<decltype(args)>(args)...));
  };
}


namespace Impl {

  // Wrapper to capture values in a lambda for perfect forwarding.
  // This captures value types by value and reference types by reference.
  template <typename T>
  struct ForwardCaptureWrapper;

  template <typename T>
  struct ForwardCaptureWrapper
  {
    template <typename TT>
    ForwardCaptureWrapper(TT&& t) : t_{std::forward<TT>(t)} {}

    auto forward() const { return std::move(t_); }

    T t_;
  };

  template <typename T>
  struct ForwardCaptureWrapper<T&>
  {
    ForwardCaptureWrapper(T& t) : t_{t} {}

    T& forward() const { return t_; };

    T& t_;
  };

  template <typename T>
  struct ForwardCaptureWrapper<const T&>
  {
    ForwardCaptureWrapper(const T& t) : t_{t} {}

    const T& forward() const { return t_; };

    const T& t_;
  };

} // end namespace Dune::Functions::Impl



/**
 * \brief Create a capture object for perfect forwarding.
 *
 * The returned object will capture the passed argument t.
 * If t is passed as r-value, then it is captured by value,
 * otherwise by reference. The captured value is accessible
 * once using the forward() method which either returns the
 * catured reference or moves the captured value.
 *
 * This allows to capture values for perfect forwarding
 * in lambda functions using
 * [t=forwardCapture(std::forward<T>(t))]() -> decltype(auto) { return t.forward(); }
 */
template <class T>
auto forwardCapture(T&& t)
{
    return Impl::ForwardCaptureWrapper<T>(std::forward<T>(t));
}



} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_COMMON_UTILITY_HH
