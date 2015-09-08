// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TUPLETREEPATH_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TUPLETREEPATH_HH


#include <tuple>


namespace Dune {
namespace Functions {


//
// Some helper functionality for tuples
namespace Std {

  // A static sequence of integers
  // just like std::integer_sequence in c++14
  template<class T, T... i>
  struct integer_sequence
  {
    using value_type = T;

    static constexpr std::size_t size() noexcept
    {
      return sizeof...(i);
    }
  };

  template<std::size_t... i>
  using index_sequence = integer_sequence<std::size_t, i...>;


} // namespace Std

namespace Imp {


  // Helper struct for construction of an
  // interger_sequence containing a range
  template<class T, T begin, T end, class Head=Std::integer_sequence<T> >
  struct MakeIntegerRangeHelper
  {};

  template<class T, T begin, T end, T... I>
  struct MakeIntegerRangeHelper<T, begin, end, Std::integer_sequence<T,I...> >
  {
    using type = typename MakeIntegerRangeHelper<T, begin+1, end, Std::integer_sequence<T, I..., begin>>::type;
  };

  template<class T, T end, T... I>
  struct MakeIntegerRangeHelper<T, end, end, Std::integer_sequence<T, I...> >
  {
    using type = Std::integer_sequence<T, I...>;
  };


  // Construct interger_sequence containg the entries begin...(end-1)
  template<class T, T begin, T end>
  using MakeIntegerRange = typename MakeIntegerRangeHelper<T, begin, end>::type;

  // Construct index_sequence containg the entries begin...(end-1)
  template<std::size_t begin, std::size_t end>
  using MakeIndexRange = MakeIntegerRange<std::size_t, begin, end>;
}



// build sub-tuple containing entries given indices in an index_sequence
template<class... T, std::size_t... i>
auto subTuple(const std::tuple<T...>& t, Std::index_sequence<i...>)
  -> decltype(std::make_tuple(std::get<i>(t)...))
{
  return std::make_tuple(std::get<i>(t)...);
}

// build sub-tuple containing a the slice consisting
// of the begin-th to (end-1)-th entries of tuple
template<std::size_t begin, std::size_t end, class... T>
auto tupleSlice(const std::tuple<T...>& t)
  -> decltype(subTuple(t, Imp::MakeIndexRange<begin, end>()))
{
  return subTuple(t, Imp::MakeIndexRange<begin, end>());
}

// build sub-tuple containing a the slice consisting
// of the first to (end-1)-th entries of tuple
template<std::size_t end, class... T>
auto tupleHead(const std::tuple<T...>& t)
  -> decltype(tupleSlice<0, end>(t))
{
  return tupleSlice<0, end>(t);
}

// build sub-tuple containing a the slice consisting
// of the begin-th to last entries of tuple
template<std::size_t begin, class... T>
auto tupleTail(const std::tuple<T...>& t)
  -> decltype(tupleSlice<begin, sizeof...(T)>(t))
{
  return tupleSlice<begin, sizeof...(T)>(t);
}





// Extend TreePath by static index given by integer i
template<int i, class... H>
std::tuple<H..., std::integral_constant<int, i>> extendTreePath(const std::tuple<H...>&head)
{
  return std::tuple_cat(head, std::make_tuple(std::integral_constant<int, i>()));
}

// Extend TreePath by dynamic index given by integer i
template<class... H>
std::tuple<H..., int> extendTreePath(const std::tuple<H...>&head, int i)
{
  return std::tuple_cat(head, std::make_tuple(i));
}





// Get child of given tree specified by given TreePath
template<class Tree>
Tree&& getChild(Tree&& tree, const std::tuple<>& tp)
{
  return tree;
}

// Get child of given tree specified by given TreePath
template<class Tree, class... I>
auto getChild(Tree&& tree, const std::tuple<int, I...>& tp)
  -> decltype(getChild(tree.child(std::get<0>(tp)), tupleTail<1>(tp)))
{
  return getChild(tree.child(std::get<0>(tp)), tupleTail<1>(tp));
}

// Get child of given tree specified by given TreePath
template<class Tree, int i0, class... I>
auto getChild(Tree&& tree, const std::tuple<std::integral_constant<int,i0>, I...>& tp)
  -> decltype(getChild(tree.template child<i0>(), tupleTail<1>(tp)))
{
  return getChild(tree.template child<i0>(), tupleTail<1>(tp));
}


} // namespace Functions
} // namespace Dune





#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TUPLETREEPATH_HH
