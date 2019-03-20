// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_INDEX_ACCESS_HH
#define DUNE_FUNCTIONS_COMMON_INDEX_ACCESS_HH


#include <utility>
#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/common/concept.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/functions/common/utility.hh>



namespace Dune {
namespace Functions {


namespace Imp {

namespace Concept {

template<class size_type>
struct HasDynamicIndexAccess
{
  template<class C>
  auto require(C&& c) -> decltype(
    c[std::declval<size_type>()]
  );
};

struct HasStaticIndexAccess
{
  template<class C>
  auto require(C&& c) -> decltype(
    c[Dune::Indices::_0]
  );
};

} // namespace Concept

} // namespace Imp



/**
 * \brief Provide operator[] index-access for containers
 *
 * \ingroup Utility
 *
 * This is the overload for types providing a operator[]
 * for dynamic std::size_t arguments.
 *
 * \param c Container to access
 * \param i The index to use for accessing the container
 * \param f A functor to call with the result of operator[]
 */
template<class C, class I, class F,
  typename std::enable_if< Dune::models<Imp::Concept::HasDynamicIndexAccess<I>, C>(), int>::type = 0>
auto hybridIndexAccess(C&& c, const I& i, F&& f)
  -> decltype(f(c[i]))
{
  return f(c[i]);
}

/**
 * \brief Provide operator[] index-access for containers
 *
 * \ingroup Utility
 *
 * This is the overload for types providing a operator[]
 * only for static arguments of type std::integral_constant<std::size_t,k>.
 * This does a static linear search until a static index
 * matching the given dynamic index is found.
 * Since the result type will in general be different
 * for different indices the method does not return
 * the result directly but passes it to a given functor.
 *
 * \param c Container to access
 * \param i The index to use for accessing the container
 * \param f A functor to call with the result of operator[]
 */
template<class C, class I, class F,
  typename std::enable_if< not Dune::models<Imp::Concept::HasDynamicIndexAccess<I>, C>(), int>::type = 0>
decltype(auto) hybridIndexAccess(C&& c, const I& i, F&& f)
{
  using Size = decltype(Hybrid::size(c));
  return Hybrid::switchCases(std::make_index_sequence<Size::value>(), i,
      [&](const auto& ii) -> decltype(auto){
        return f(c[ii]);
      }, [&]() -> decltype(auto){
        return f(c[Dune::Indices::_0]);
      });
}


namespace Imp {

  /**
   * \brief Class representing a shifted multi index
   *
   * \tparam Index Type of the base multi index
   * \tparam offset Number of positions to shift left
   *
   * For a given multi index of size n this
   * represents a multi index with the first
   * offset entries removed.
   *
   * Notice that this does only store a reference to
   * the passed multi index.
   */
  template<class Index, std::size_t offset=1>
  class ShiftedDynamicMultiIndex
  {
  public:
    ShiftedDynamicMultiIndex(const Index& index) :
      index_(index)
    {}

    std::size_t operator[](std::size_t position) const
    {
      if (position<size())
        return index_[position+offset];
      else
        return 0;
    }

    /**
     * \brief Return multi index with one more position truncated
     */
    ShiftedDynamicMultiIndex<Index, offset+1> pop() const
    {
      return {index_};
    }

    std::size_t size() const
    {
      if (offset < index_.size())
        return index_.size() - offset;
      else
        return 0;
    }

  private:
    const Index& index_;
  };

  template<class Index, std::size_t offset=1>
  class ShiftedStaticMultiIndex
  {
  public:
    ShiftedStaticMultiIndex(const Index& index) :
      index_(index)
    {}

    template<std::size_t i>
    auto operator[](Dune::index_constant<i>) const
    {
      auto isContained = Dune::Std::bool_constant<(i<size())>{};
      return Hybrid::ifElse(isContained, [&](auto id) {
        return id(index_)[Dune::index_constant<i+offset>{}];
      }, [](auto id) {
        return Dune::Indices::_0;
      });
    }

    /**
     * \brief Return multi index with one more position truncated
     */
    ShiftedStaticMultiIndex<Index, offset+1> pop() const
    {
      return {index_};
    }

    static constexpr std::size_t size()
    {
      auto fullSize = decltype(Hybrid::size(std::declval<Index>()))::value;
      if (offset < fullSize)
        return fullSize - offset;
      else
        return 0;
    }

  private:
    const Index& index_;
  };

  /**
   * \brief Create a ShiftedDynamicMultiIndex
   *
   * \tparam offset Number of positions to shift left
   */
  template<std::size_t offset, class Index>
  ShiftedDynamicMultiIndex<Index, offset> shiftedDynamicMultiIndex(const Index& index)
  {
    return {index};
  }

  template<std::size_t offset, class Index>
  ShiftedStaticMultiIndex<Index, offset> shiftedStaticMultiIndex(const Index& index)
  {
    return {index};
  }

} // namespace Imp




namespace Imp {

template<class Result, class Index>
struct MultiIndexResolver
{
  MultiIndexResolver(const Index& index) :
    index_(index)
  {}

  template<class C,
    typename std::enable_if<not std::is_convertible<C&, Result>::value, int>::type = 0>
  Result operator()(C&& c)
  {
    auto&& subIndex = Imp::shiftedDynamicMultiIndex<1>(index_);
    auto&& subIndexResolver = MultiIndexResolver<Result, decltype(subIndex)>(subIndex);
    return (Result)(hybridIndexAccess(c, index_[Dune::Indices::_0], subIndexResolver));
  }

  template<class C,
    typename std::enable_if<std::is_convertible<C&, Result>::value, int>::type = 0>
  Result operator()(C&& c)
  {
    return (Result)(std::forward<C>(c));
  }

  const Index& index_;
};

} // namespace Imp



/**
 * \brief Provide multi-index access by chaining operator[]
 *
 * \ingroup Utility
 *
 * This provides access to a nested container by given
 * multi-index. Internally this is resolved by recusive
 * operator[]-calls with static or dynamic indices.
 * Because this recursion must be terminated using a
 * compile-time criterion, the result type must explicitly
 * be provided. The recursion will terminate once the
 * result can be converted to this result type.
 *
 * \tparam Result Type of result
 *
 * \param c Container to access
 * \param index Multi-index
 */
template<class Result, class C, class MultiIndex>
Result hybridMultiIndexAccess(C&& c, const MultiIndex& index)
{

  Imp::MultiIndexResolver<Result, MultiIndex> multiIndexResolver(index);
  return multiIndexResolver(c);
}






namespace Imp {

  template<class C, class MultiIndex, class IsFinal>
  constexpr decltype(auto) resolveDynamicMultiIndex(C&& c, const MultiIndex& multiIndex, const IsFinal& isFinal)
  {
    // If c is already considered final simply return it,
    // else resolve the next multiIndex entry.
    return Hybrid::ifElse(isFinal(c), [&, c = forwardCapture(std::forward<C>(c))](auto id) -> decltype(auto) {
      assert(multiIndex.size() == 0);
      return c.forward();
    }, [&](auto id) -> decltype(auto) {
      auto hasDynamicAccess = callableCheck([](auto&& cc) -> void_t<decltype(cc[0])> {});

      // Split multiIndex into first entry and remaining ones.
      auto i = multiIndex[0];
      auto tail = multiIndex.pop();

      // Resolve first multiIndex entry by c[multiIndex[0]] and
      // continue resolving with the remaining remaining ones.
      // If c has a dynamic operator[] this is straight forward.
      // Else the dynamic multiIndex[0] has to be translated into
      // a static one using hybridIndexAccess.
      return Hybrid::ifElse(hasDynamicAccess(c), [&](auto id) -> decltype(auto) {
        return Imp::resolveDynamicMultiIndex(id(c)[i], tail, isFinal);
      }, [&](auto id) -> decltype(auto) {
        // auto indexRange = range(Hybrid::size(id(c)));
        auto indexRange = typename decltype(range(Hybrid::size(id(c))))::integer_sequence();
        return Hybrid::switchCases(indexRange, i, [&](auto static_i) -> decltype(auto){
          // Do rescursion with static version of i
          return Imp::resolveDynamicMultiIndex(id(c)[static_i], tail, isFinal);
        }, [&]() -> decltype(auto){
          // As fallback we use c[0] this is needed, because there must be one branch that matches.
          return Imp::resolveDynamicMultiIndex(id(c)[Dune::Indices::_0], tail, isFinal);
        });
      });
    });
  }

  template<class C, class MultiIndex>
  constexpr decltype(auto) resolveStaticMultiIndex(C&& c, const MultiIndex& multiIndex)
  {
    auto isExhausted = Hybrid::equals(Hybrid::size(multiIndex), Dune::Indices::_0);
    return Hybrid::ifElse(isExhausted, [&, c = forwardCapture(std::forward<C>(c))](auto id) -> decltype(auto) {
      return c.forward();
    }, [&](auto id) -> decltype(auto) {
      auto head = multiIndex[Dune::Indices::_0];
      auto tail = multiIndex.pop();

      return Imp::resolveStaticMultiIndex(id(c)[head], tail);
    });
  }

} // namespace Imp



/**
 * \brief Provide multi-index access by chaining operator[]
 *
 * \ingroup Utility
 *
 * This provides access to a nested container by given dynamically sized
 * multi-index. Internally this is resolved by recusive
 * operator[]-calls with static or dynamic indices.
 * Because this recursion must be terminated, a predicate
 * is used to determine if a type is considered final for the
 * multi-index resolution. Hence multi-index resolution is
 * terminated for values where the criterion matches.
 * In order to do this statically the predicate has to
 * return its result as std::true_type or std::false_type.
 *
 * If the entries of the multi-index are exhausted, additional
 * [0] entries are used until the termination criterion is satisfied.
 *
 * \param c Container to access
 * \param index Multi-index
 * \param isFinal A predicate function checking if recursion should be terminated.
 */
template<class C, class MultiIndex, class IsFinal>
constexpr decltype(auto) resolveDynamicMultiIndex(C&& c, const MultiIndex& multiIndex, const IsFinal& isFinal)
{
  return Imp::resolveDynamicMultiIndex(std::forward<C>(c), Imp::shiftedDynamicMultiIndex<0>(multiIndex), isFinal);
}

/**
 * \brief Provide multi-index access by chaining operator[]
 *
 * \ingroup Utility
 *
 * This provides access to a nested container by given dynamically sized
 * multi-index. Internally this is resolved by recusive
 * operator[]-calls with static or dynamic indices.
 * Because this recursion must be terminated, a predicate
 * is used to determine if a type is considered final for the
 * multi-index resolution. This version terminates the recursion
 * on values that neither have a static nore a dynamic operator[].
 *
 * \param c Container to access
 * \param index Multi-index
 */
template<class C, class MultiIndex>
constexpr decltype(auto) resolveDynamicMultiIndex(C&& c, const MultiIndex& multiIndex)
{
  auto hasNoIndexAccess = negatePredicate(callableCheck([](auto&& cc) -> void_t<decltype(cc[Dune::Indices::_0])> {}));
  return Imp::resolveDynamicMultiIndex(std::forward<C>(c), Imp::shiftedDynamicMultiIndex<0>(multiIndex), hasNoIndexAccess);
}

/**
 * \brief Provide multi-index access by chaining operator[]
 *
 * \ingroup Utility
 *
 * This provides access to a nested container by given statically sized
 * multi-index. Internally this is resolved by recusive
 * operator[]-calls with static or dynamic indices.
 * Since the multi-index must have compile-time known size
 * it is possible, that values resolved by different multi-indices
 * have a different size.
 *
 * \param c Container to access
 * \param index Multi-index
 */
template<class C, class MultiIndex>
constexpr decltype(auto) resolveStaticMultiIndex(C&& c, const MultiIndex& multiIndex)
{
  return Imp::resolveStaticMultiIndex(std::forward<C>(c), Imp::shiftedStaticMultiIndex<0>(multiIndex));
}



} // namespace Dune::Functions
} // namespace Dune



#endif // DUNE_FUNCTIONS_COMMON_INDEX_ACCESS_HH
