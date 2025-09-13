// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH

#include <array>
#include <cassert>
#include <functional>
#include <type_traits>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/tuplevector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>

/**
 * \file containerdescriptor.hh
 * \brief Lightweight representation of (hierarchical) size and block structure extracted
 * from a basis to describe data structures like containers that can be accessed by
 * multi-indices provided by the basis.
 *
 * The structure of a container-descriptor is a reduced container interface:
 * \code
  struct [Container]Descriptor
  {
    template<class Index>
    [SubContainerDescriptor] operator[](Index i) const;  // return the i-th sub-container-descriptor

    [static constexpr] std::size_t size() [const];  // return the number of children
  };
 * \endcode
 *
 * With the `operator[]` you can access the children. The `Index` type is either
 * an integral value or an `integral_constant` for tuple nodes.
 *
 * Size is either a static property, or a runtime value.
 **/

namespace Dune::Functions {
namespace ContainerDescriptors {

//! Fallback container descriptor if nothing else fits
struct Unknown {};

} // end namespace ContainerDescriptors

//! Return the container descriptor of the pre-basis, if defined, otherwise ContainerDescriptor::Unknown
template<class PreBasis>
auto containerDescriptor(const PreBasis& preBasis)
{
  if constexpr (requires{ preBasis.containerDescriptor(); })
    return preBasis.containerDescriptor();
  else
    return ContainerDescriptors::Unknown{};
}


namespace ContainerDescriptors {

//! The node in the descriptor tree representing a value placeholder
struct Value
{
  //! The child access method is only available for the interface, but should not be called.
  template<class Index>
  Value operator[] (const Index&) const { return {}; }

  //! A value placeholder does not have any sub-descriptors, thus its size is zero.
  static constexpr std::size_t size () { return 0; }
};

//! Descriptor with all children of possibly different type
template<class... Children>
using Tuple = Dune::TupleVector<Children...>;

//! Generate a descriptor in case the children are not all of the same type.
//! \relates Tuple
template<class Child0, class... Children,
  std::enable_if_t<(sizeof...(Children) > 0), int> = 0,
  std::enable_if_t<(...|| (not std::is_same_v<Child0, Children>)), int> = 0>
auto makeDescriptor (Child0 child0, Children... children)
{
  using Descriptor = Tuple<Child0,Children...>;
  return Descriptor{std::move(child0),std::move(children)...};
}


//! Descriptor for arrays with all children of the same type and static size.
template<class Child, std::size_t n>
using Array = std::array<Child, n>;

//! Generate a descriptor in case the children are all of the same type.
template<class Child0, class... Children,
  std::enable_if_t<(std::is_same_v<Child0, Children> &&...), int> = 0>
auto makeDescriptor (Child0 child, Children... children)
{
  using Descriptor = Array<Child0,1+sizeof...(Children)>;
  return Descriptor{std::move(child),std::move(children)...};
}


//! Descriptor for vectors with all children of the same type and dynamic size.
template<class Child>
using Vector = std::vector<Child>;

//! Descriptor for arrays with all children identical and the number of children a static size.
template<class Child, std::size_t n>
struct UniformArray
{
  //! Default constructor. Is enable if the child-type is default constructible.
  template<class C = Child,
    std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
  UniformArray ()
    : child_{}
  {}

  //! Constructor that stores a single child only.
  explicit UniformArray (Child child)
    : child_{std::move(child)}
  {}

  //! Access the i'th child that is always the same, i.e., `child_`.
  template<class Index>
  const Child& operator[] (const Index& /*i*/) const { return child_; }

  //! The static size information, i.e., number of children.
  static constexpr std::size_t size () { return n; }

private:
  Child child_;
};

//! Alias for a uniform array storing value placeholders
template<std::size_t n>
using FlatArray = UniformArray<Value,n>;

//! Generate a uniform descriptor in case the size is a static constant
template<class Child, std::size_t n>
auto makeUniformDescriptor (std::integral_constant<std::size_t,n>, Child child)
{
  return UniformArray<Child,n>{std::move(child)};
}


//! Uniform descriptor with dynamic size.
template<class Child>
struct UniformVector
{
  //! Default constructor with size. Is enable if the child-type is default constructible.
  template<class C = Child,
    std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
  explicit UniformVector (std::size_t size)
    : size_{size}
    , child_{}
  {}

  //! Constructor that stores the size and a single child only.
  UniformVector (std::size_t size, Child child)
    : size_{size}
    , child_{std::move(child)}
  {}

  //! Access the i'th child that is always the same, i.e., `child_`.
  template<class Index>
  const Child& operator[] (const Index& /*i*/) const { return child_; }

  //! The dynamic size information, i.e., number of children.
  std::size_t size () const { return size_; }

private:
  std::size_t size_;
  Child child_;
};

//! Alias for a uniform vector storing value placeholders
using FlatVector = UniformVector<Value>;

//! Generate a uniform descriptor in case the size is a dynamic value
template<class Child>
auto makeUniformDescriptor (std::size_t n, Child child)
{
  return UniformVector<Child>{n,std::move(child)};
}

namespace Impl {

template<class InnerFunc, class LeafFunc>
struct TreeTransform
{
  TreeTransform (const InnerFunc& innerFunc, const LeafFunc& leafFunc)
    : innerFunc_(innerFunc)
    , leafFunc_(leafFunc)
  {}

  Unknown operator() (const Unknown& tree) const
  {
    return tree;
  }

  auto operator() (const Value& tree) const
  {
    return leafFunc_(tree);
  }

  template<class... V>
  auto operator() (const Tuple<V...>& tree) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return makeDescriptor(innerFunc_(tree[ii])...);
    }, std::make_index_sequence<sizeof...(V)>());
  }

  template<class V, std::size_t n>
  auto operator() (const Array<V,n>& tree) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return makeDescriptor(innerFunc_(tree[ii])...);
    }, std::make_index_sequence<n>());
  }

  template<class V>
  auto operator() (const Vector<V>& tree) const
  {
    using W = decltype(innerFunc_(tree[0]));
    Vector<W> result;
    result.reserve(tree.size());
    for (std::size_t i = 0; i < tree.size(); ++i)
      result.emplace_back(innerFunc_(tree[i]));
    return result;
  }

  template<class V, std::size_t n>
  auto operator() (const UniformArray<V,n>& tree) const
  {
    return makeUniformDescriptor(Dune::index_constant<n>{}, innerFunc_(tree[0]));
  }

  template<class V>
  auto operator() (const UniformVector<V>& tree) const
  {
    return makeUniformDescriptor(tree.size(), innerFunc_(tree[0]));
  }

private:
  InnerFunc innerFunc_;
  LeafFunc leafFunc_;
};


/**
  * Append a size to the leaf-nodes of the tree
  *
  * This transform of the given tree is used to implement
  * a blocked-interleaved index-merging strategy in a power-basis.
  *
  * Examples:
  * append( Flat[Container] c, size ) -> Uniform[Container]( c.size(), Flat[Container](size) )
  * append( Descriptor(child...), size ) -> Descriptor( append(child, size)... )
  */
template<class Size, class T>
auto appendToTree (Size s, const T& tree)
{
  auto transform = TreeTransform(
    [s](auto&& node) { return appendToTree(s, node); },
    [s](auto&& node) { return makeUniformDescriptor(s, node); });
  return transform(tree);
}

//! Concatenate n descriptors by flat lexicographic merging wrt first digit
template<class... Child>
auto flatLexicographic (Child... child)
{
  if constexpr ((std::is_same_v<Child, FlatVector> && ...))
    return FlatVector((child.size() + ...));
  else
    return Unknown{};
}

//! Concatenate descriptor n-times by flat lexicographic merging wrt first digit
template<class N, class Child>
auto flatLexicographicN (N n, Child child)
{
  return Unknown{};
}

//! Concatenate descriptor n-times by flat lexicographic merging wrt first digit
template<class N, class GrandChild>
auto flatLexicographicN (N n, UniformVector<GrandChild> child)
{
  return UniformVector<GrandChild>{child.size()*n, child[0]};
}

//! Concatenate descriptor n-times by flat lexicographic merging wrt first digit
template<class N, class GrandChild, std::size_t m>
auto flatLexicographicN (N n, UniformArray<GrandChild, m> child)
{
  if constexpr (std::is_same_v<N, std::size_t>)
    return UniformVector<GrandChild>{n*m, child[0]};
  else
    return UniformArray<GrandChild, N::value*m>{child[0]};
}

//! Concatenate descriptor n-times by flat lexicographic merging wrt first digit
template<class N, class GrandChild>
auto flatLexicographicN (N n, Vector<GrandChild> child)
{
  auto result = Vector<GrandChild>{};
  result.reserve(child.size()*n);
  for (auto j : Dune::range(n))
    for (auto i : Dune::range(child.size()))
      result.emplace_back(child[i]);
  return result;
}

//! Concatenate descriptor n-times by flat lexicographic merging wrt first digit
template<class N, class GrandChild, std::size_t m>
auto flatLexicographicN (N n, Array<GrandChild, m> child)
{
  if constexpr (std::is_same_v<N, std::size_t>)
  {
    auto result = Vector<GrandChild>{};
    result.reserve(child.size()*n);
    for (auto j : Dune::range(n))
      for (auto i : Dune::range(child.size()))
        result.emplace_back(child[i]);
    return result;
  }
  else
  {
    auto result = Array<GrandChild, N::value*m>{};
    for (auto j : Dune::range(n))
      for (auto i : Dune::range(child.size()))
        result.emplace_back(child[i]);
    return result;
  }
}

//! Concatenate descriptor n-times by flat lexicographic merging wrt first digit
template<class N, class... GrandChild>
auto flatLexicographicN (N n, Tuple<GrandChild...> child)
{
  constexpr std::size_t m = sizeof...(GrandChild);
  return Dune::unpackIntegerSequence([&](auto... i) {
    return makeDescriptor(std::get<i%m>(child)...);
  }, std::make_index_sequence<n*m>{});
}

//! Concatenate descriptor n-times by flat interleaved merging wrt first digit
template<class N, class Child>
auto flatInterleavedN (N n, Child child)
{
  return Unknown{};
}

//! Concatenate descriptor n-times by flat interleaved merging wrt first digit
template<class N, class GrandChild>
auto flatInterleavedN (N n, UniformVector<GrandChild> child)
{
  return UniformVector<GrandChild>{child.size()*n, child[0]};
}

//! Concatenate descriptor n-times by flat interleaved merging wrt first digit
template<class N, class GrandChild, std::size_t m>
auto flatInterleavedN (N n, UniformArray<GrandChild, m> child)
{
  if constexpr (std::is_integral_v<N>)
    return UniformVector<GrandChild>{n*m, child[0]};
  else
    return UniformArray<GrandChild, N::value*m>{child[0]};
}

//! Concatenate descriptor n-times by flat interleaved merging wrt first digit
template<class N, class GrandChild>
auto flatInterleavedN (N n, Vector<GrandChild> child)
{
  auto result = Vector<GrandChild>{};
  result.reserve(child.size()*n);
  for (auto i : Dune::range(child.size()))
    for (auto j : Dune::range(n))
      result.emplace_back(child[i]);
  return result;
}

//! Concatenate descriptor n-times by flat interleaved merging wrt first digit
template<class N, class GrandChild, std::size_t m>
auto flatInterleavedN (N n, Array<GrandChild, m> child)
{
  if constexpr (std::is_integral_v<N>)
  {
    auto result = Vector<GrandChild>{};
    result.reserve(child.size()*n);
    for (auto i : Dune::range(child.size()))
      for (auto j : Dune::range(n))
        result.emplace_back(child[i]);
    return result;
  }
  else
  {
    auto result = Array<GrandChild, N::value*m>{};
    for (auto i : Dune::range(child.size()))
      for (auto j : Dune::range(n))
        result.emplace_back(child[i]);
    return result;
  }
}

//! Concatenate descriptor n-times by flat interleaved merging wrt first digit
template<class N, class... GrandChild>
auto flatInterleavedN (N n, Tuple<GrandChild...> child)
{
  constexpr std::size_t m = sizeof...(GrandChild);
  return Dune::unpackIntegerSequence([&](auto... i) {
    return makeDescriptor(std::get<i/n>(child)...);
  }, std::make_index_sequence<n*m>{});
}

} // end namespace Impl
} // end namespace ContainerDescriptors
} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH
