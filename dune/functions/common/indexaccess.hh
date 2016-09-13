// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_INDEX_ACCESS_HH
#define DUNE_FUNCTIONS_COMMON_INDEX_ACCESS_HH


#include <dune/common/concept.hh>
#include <dune/common/hybridutilities.hh>



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
class ShiftedMultiIndex
{
public:
  ShiftedMultiIndex(const Index& index) :
    index_(index)
  {}

  template<class P>
  decltype(auto) operator[](const P& position) const
  {
    return index_[position+offset];
  }

  /**
   * \brief Return multi index with one more position truncated
   */
  ShiftedMultiIndex<Index, offset+1> pop() const
  {
    return {index_};
  }

  auto size() const
  {
    return index_.size() - offset;
  }

private:
  const Index& index_;
};



/**
 * \brief Create a ShiftedMultiIndex
 *
 * \tparam offset Number of positions to shift left
 */
template<std::size_t offset, class Index>
ShiftedMultiIndex<Index, offset> shiftedMultiIndex(const Index& index)
{
  return {index};
}



/**
 * \brief Create a ShiftedMultiIndex with one position truncated
 *
 * \tparam offset Number of positions to shift left
 */
template<class Index>
ShiftedMultiIndex<Index, 1> shiftedMultiIndex(const Index& index)
{
  return {index};
}



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
    auto&& subIndex = shiftedMultiIndex(index_);
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



} // namespace Dune::Functions
} // namespace Dune



#endif // DUNE_FUNCTIONS_COMMON_INDEX_ACCESS_HH
