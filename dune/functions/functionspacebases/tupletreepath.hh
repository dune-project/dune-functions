// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TUPLETREEPATH_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TUPLETREEPATH_HH


#include <tuple>
#include <dune/common/std/utility.hh>


namespace Dune {
namespace Functions {

namespace StaticIndices {
  static const std::integral_constant<int, 0> _0 =
    std::integral_constant<int, 0>();
  static const std::integral_constant<int, 1> _1 =
    std::integral_constant<int, 1>();
  static const std::integral_constant<int, 2> _2 =
    std::integral_constant<int, 2>();
  static const std::integral_constant<int, 3> _3 =
    std::integral_constant<int, 3>();
  static const std::integral_constant<int, 4> _4 =
    std::integral_constant<int, 4>();
  static const std::integral_constant<int, 5> _5 =
    std::integral_constant<int, 5>();
  static const std::integral_constant<int, 6> _6 =
    std::integral_constant<int, 6>();
  static const std::integral_constant<int, 7> _7 =
    std::integral_constant<int, 7>();
  static const std::integral_constant<int, 8> _8 =
    std::integral_constant<int, 8>();
  static const std::integral_constant<int, 9> _9 =
    std::integral_constant<int, 9>();
}


//
// Some helper functionality for tuples
namespace Std {

  // A static sequence of indices
  // just like std::index_sequence in c++14
  //
  // While we already have Dune::Std::index_sequence<...>, there's a
  // major problem with it. The latter does only derive form
  // Dune::Std::integer_sequence<std::size_t,...> but is not
  // equal to. Hence all template specializations for integer_sequence
  // (as done e.g. in the utilities below) would not work without
  // special treatment for index_sequence. The proper implementation
  // is given by the following alias in Dune::Functions::Std.
  template<std::size_t... i>
  using index_sequence = Dune::Std::integer_sequence<std::size_t, i...>;

}


namespace Imp
{

  // Helper struct for construction of an
  // interger_sequence containing a range
  template<class T, T begin, T end, class Head=Dune::Std::integer_sequence<T> >
  struct MakeIntegerRangeHelper
  {};

  template<class T, T begin, T end, T... I>
  struct MakeIntegerRangeHelper<T, begin, end, Dune::Std::integer_sequence<T,I...> >
  {
    using type = typename MakeIntegerRangeHelper<T, begin+1, end, Dune::Std::integer_sequence<T, I..., begin>>::type;
  };

  template<class T, T end, T... I>
  struct MakeIntegerRangeHelper<T, end, end, Dune::Std::integer_sequence<T, I...> >
  {
    using type = Dune::Std::integer_sequence<T, I...>;
  };

}

// Construct interger_sequence containg the entries begin...(end-1)
template<class T, T begin, T end>
using MakeIntegerRange = typename Imp::MakeIntegerRangeHelper<T, begin, end>::type;

// Construct index_sequence containg the entries begin...(end-1)
template<std::size_t begin, std::size_t end>
using MakeIndexRange = MakeIntegerRange<std::size_t, begin, end>;



// build sub-tuple containing entries given indices in an index_sequence
template<class... T, std::size_t... i>
auto subTuple(const std::tuple<T...>& t, Functions::Std::index_sequence<i...>)
  -> decltype(std::make_tuple(std::get<i>(t)...))
{
  return std::make_tuple(std::get<i>(t)...);
}

// build sub-tuple containing a the slice consisting
// of the begin-th to (end-1)-th entries of tuple
template<std::size_t begin, std::size_t end, class... T>
auto tupleSlice(const std::tuple<T...>& t)
  -> decltype(subTuple(t, MakeIndexRange<begin, end>()))
{
  return subTuple(t, MakeIndexRange<begin, end>());
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





int canonicalizeTreePathIndex(int i)
{
  return i;
}

template<int i>
std::integral_constant<int,i> canonicalizeTreePathIndex(const std::integral_constant<int, i>& ci)
{
  return std::integral_constant<int,i>();
}

// Create a TreePath from arguments
template<class... I>
auto makeTreePath(I&&... i)
  -> decltype(std::make_tuple(canonicalizeTreePathIndex(i)...))
{
  return std::make_tuple(canonicalizeTreePathIndex(i)...);
}

// Extend TreePath by static index given by integer i
template<int i, class... H>
std::tuple<H..., std::integral_constant<int, i>> extendTreePath(const std::tuple<H...>&head, typename std::integral_constant<int,i> I ={})
{
  return std::tuple_cat(head, std::make_tuple(I));
}

// Extend TreePath by dynamic index given by integer i
template<class... H>
std::tuple<H..., int> extendTreePath(const std::tuple<H...>&head, int i)
{
  return std::tuple_cat(head, std::make_tuple(i));
}



// Get child with empty TreePath, i.e., the tree itself
template<class Tree>
auto getChild(Tree&& tree)
  ->decltype(std::forward<Tree>(tree))
{
  return std::forward<Tree>(tree);
}

// Get child of given tree specified by given TreePath
template<class Tree, class... I>
auto getChild(Tree&& tree, const int& i0, const I&... i)
  -> decltype(getChild(tree.child(i0), i...))
{
  return getChild(tree.child(i0), i...);
}

template<class Tree, int i0, class... I>
auto getChild(Tree&& tree, const std::integral_constant<int,i0>&, const I&... i)
  -> decltype(getChild(tree. template child<i0>(), i...))
{
  return getChild(tree. template child<i0>(), i...);
}

namespace Imp {

  template<class Tree, class... T, std::size_t... i>
  auto getChildHelper(Tree&& tree, const std::tuple<T...>& tp, const Functions::Std::index_sequence<i...>&)
    -> decltype(getChild(tree, std::get<i>(tp)...))
  {
    return getChild(tree, std::get<i>(tp)...);
  }

}

// Get child of given tree specified by given TreePath
template<class Tree, class... I>
auto getChild(Tree&& tree, const std::tuple<I...>& tp)
  -> decltype(Imp::getChildHelper(tree, tp, MakeIndexRange<0,sizeof...(I)>()))
{
  return Imp::getChildHelper(tree, tp, MakeIndexRange<0,sizeof...(I)>());
}



namespace Imp {

template<class Node, int... Path>
struct ChildTypeHelper{};

template<class Node, int FirstChild, int... Path>
struct ChildTypeHelper<Node, FirstChild, Path...>
{
  using Type = typename ChildTypeHelper<typename Node::template Child<FirstChild>::Type, Path...>::Type;
};

template<class Node>
struct ChildTypeHelper<Node>
{
  using Type = Node;
};


}

// Simplify access to type of child
template<class Node, int... Path>
using ChildType = typename Imp::ChildTypeHelper<Node, Path...>::Type;


} // namespace Functions
} // namespace Dune





#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TUPLETREEPATH_HH
