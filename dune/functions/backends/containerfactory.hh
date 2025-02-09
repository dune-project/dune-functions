// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_BACKENDS_CONTAINERFACTORY_HH
#define DUNE_FUNCTIONS_BACKENDS_CONTAINERFACTORY_HH

#include <type_traits>
#include <array>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/tuplevector.hh>

#include <dune/functions/functionspacebases/containerdescriptors.hh>

namespace Dune::Functions {
namespace ContainerDescriptors {
namespace Impl {

template<class T>
struct ContainerFactory
{
  void operator() (const Unknown& descriptor, const T& defaultValue) const
  {
    DUNE_THROW(Dune::NotImplemented, "Cannot create a vector. The container descriptor is unknown.");
  }

  template<class... ChildDescriptor>
  auto operator() (const Tuple<ChildDescriptor...>& descriptor, const T& defaultValue) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return Dune::TupleVector<decltype((*this)(descriptor[ii], defaultValue))...>{(*this)(descriptor[ii], defaultValue)...};
    }, std::make_index_sequence<sizeof...(ChildDescriptor)>());
  }

  template<class ChildDescriptor, std::size_t n>
  auto operator() (const Array<ChildDescriptor,n>& descriptor, const T& defaultValue) const
  {
    using ChildContainer = decltype((*this)(descriptor[0], defaultValue));
    return unpackIntegerSequence([&](auto... ii) {
      return std::array<ChildContainer, n>{(*this)(descriptor[ii], defaultValue)...};
    }, std::make_index_sequence<n>());
  }

  template<class ChildDescriptor>
  auto operator() (const Vector<ChildDescriptor>& descriptor, const T& defaultValue) const
  {
    using ChildContainer = decltype((*this)(descriptor[0], defaultValue));
    auto result = std::vector<ChildContainer>();
    result.reserve(descriptor.size());
    for (std::size_t i = 0; i < descriptor.size(); ++i)
      result.emplace_back((*this)(descriptor[i], defaultValue));
    return result;
  }

  template<class ChildDescriptor, std::size_t n>
  auto operator() (const UniformArray<ChildDescriptor,n>& descriptor, const T& defaultValue) const
  {
    using ChildContainer = decltype((*this)(descriptor[0], defaultValue));
    auto childContainer = (*this)(descriptor[0], defaultValue);
    return unpackIntegerSequence([&](auto... ii) {
      return std::array<ChildContainer, n>{((void)(ii),childContainer)...};
    }, std::make_index_sequence<n>());
  }

  template<class ChildDescriptor>
  auto operator() (const UniformVector<ChildDescriptor>& descriptor, const T& defaultValue) const
  {
    using ChildContainer = decltype((*this)(descriptor[0], defaultValue));
    auto childContainer = (*this)(descriptor[0], defaultValue);
    return std::vector<ChildContainer>(descriptor.size(), childContainer);
  }

  // scalar types

  auto operator() (const Value& descriptor, const T& defaultValue) const
  {
    return T(defaultValue);
  }

};

} // end namespace Impl
} // end namespace ContainerDescriptors


/**
 * \brief Construct a nested random access container compatible with the container descriptor.
 *
 * \param descriptor A ContainerDescriptor provided by a global basis
 * \param defaultValue The default value to initialize the container entries.
 *
 * The constructed container mimics the nested structure of the container descriptor,
 * but uses data structures like `std::vector`, `std::array`, and `Dune::TupleVector`
 * to represent the block levels. The entries in the vector are of type `T` and
 * initialized with the provided default value.
 */
template<class T, class ContainerDescriptor>
auto makeContainer (const ContainerDescriptor& descriptor, const T& defaultValue)
{
  auto factory = ContainerDescriptors::Impl::ContainerFactory<T>{};
  return factory(descriptor, defaultValue);
}

/**
 * \brief Construct a nested random access container compatible with the container descriptor.
 *
 * \param descriptor A ContainerDescriptor provided by a global basis
 *
 * The constructed container mimics the nested structure of the container descriptor,
 * but uses data structures like `std::vector`, `std::array`, and `Dune::TupleVector`
 * to represent the block levels. The entries in the vector are of type `T` and
 * default initialized.
 */
template<class T, class ContainerDescriptor>
auto makeContainer (const ContainerDescriptor& descriptor)
{
  return makeContainer<T>(descriptor, T());
}

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_BACKENDS_CONTAINERFACTORY_HH
