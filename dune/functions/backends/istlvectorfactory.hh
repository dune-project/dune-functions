// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_BACKENDS_ISTL_VECTORFACTORY_HH
#define DUNE_FUNCTIONS_BACKENDS_ISTL_VECTORFACTORY_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>

#include <dune/functions/functionspacebases/containerdescriptors.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

namespace Dune::Functions {
namespace ContainerDescriptors {

template<class T>
struct ISTLVectorFactory
{
  void operator() (const Unknown& tree) const
  {
    DUNE_THROW(Dune::NotImplemented, "Cannot create a vector. The container descriptor is unknown.");
  }

  template<class... V>
  auto operator() (const Tuple<V...>& tree) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return Dune::MultiTypeBlockVector<decltype((*this)(tree[ii]))...>{(*this)(tree[ii])...};
    }, std::make_index_sequence<sizeof...(V)>());
  }

  template<class V, std::size_t n>
  auto operator() (const Array<V,n>& tree) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return Dune::BlockVector{(*this)(tree[ii])...};
    }, std::make_index_sequence<n>());
  }

  template<class V>
  auto operator() (const Vector<V>& tree) const
  {
    using W = decltype((*this)(tree[0]));
    Dune::BlockVector<W> result(tree.size());
    for (std::size_t i = 0; i < tree.size(); ++i)
      result[i] = (*this)(tree[i]);
    return result;
  }

  template<class V, std::size_t n>
  auto operator() (const UniformArray<V,n>& tree) const
  {
    auto node = (*this)(tree[0]);
    return unpackIntegerSequence([&](auto... ii) {
      return Dune::BlockVector{((void)(ii),node)...};
    }, std::make_index_sequence<n>());
  }

  template<class V>
  auto operator() (const UniformVector<V>& tree) const
  {
    auto node = (*this)(tree[0]);
    using W = decltype(node);
    Dune::BlockVector<W> result(tree.size());
    for (std::size_t i = 0; i < tree.size(); ++i)
      result[i] = node;
    return result;
  }

  // scalar types

  auto operator() (const Value& tree) const
  {
    return T(0);
  }

  // flat vectors

  template<std::size_t n>
  auto operator() (const UniformArray<Value,n>& tree) const
  {
    return Dune::FieldVector<T,n>(0);
  }

  auto operator() (const UniformVector<Value>& tree) const
  {
    return Dune::BlockVector<T>(tree.size());
  }

  // block vectors

  template<std::size_t n>
  auto operator() (const UniformVector<UniformArray<Value,n>>& tree) const
  {
    return Dune::BlockVector<Dune::FieldVector<T,n>>(tree.size());
  }

  template<std::size_t n>
  auto operator() (const Vector<UniformArray<Value,n>>& tree) const
  {
    return Dune::BlockVector<Dune::FieldVector<T,n>>(tree.size());
  }

  template<std::size_t n, std::size_t m>
  auto operator() (const Array<UniformArray<Value,n>,m>& tree) const
  {
    return Dune::BlockVector<Dune::FieldVector<T,n>>(m);
  }
};

} // end namespace ContainerDescriptors


// Construct an istl vector type compatible with the container descriptor
template<class T = double, class ContainerDescriptor>
auto istlVectorFactory (const ContainerDescriptor& tree)
{
  auto factory = ContainerDescriptors::ISTLVectorFactory<T>{};
  return factory(tree);
}

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_BACKENDS_ISTL_VECTORFACTORY_HH
