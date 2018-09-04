#pragma once

#include <algorithm>
#include <memory>
#include <type_traits>

namespace Dune {
namespace Functions {

// TODO IMHO this should move to dune-common ASAP

namespace Impl {

//! Wrap const reference as shared_ptr
template <class T, class S>
std::shared_ptr<const T> wrap_own_share(const S& t) {
  return std::shared_ptr<T>(&t, [](auto*) {});
}

//! Wrap reference as shared_ptr
template <class T, class S>
std::shared_ptr<T> wrap_own_share(S& t) {
  return std::shared_ptr<T>(&t, [](auto*) {});
}

//! Move r-value reference to shared_ptr
template <class T, class S>
std::shared_ptr<T> wrap_own_share(S&& t) {
  return std::make_shared<S>(std::move(t));
}

//! Share ownership of shared_ptr
template <class T, class S>
std::shared_ptr<T> wrap_own_share(std::shared_ptr<S> t) {
  return t;
}

// TODO move to generic include
template <class T>
using DisablePointer = std::enable_if_t<not std::is_pointer<T>::value>;

} // end namespace Impl

/** \brief Convert l-value and r-value references, and shared_ptr of an object
 * into a shared_ptr of a convertible type.
 *
 * L-value references are wrapped, r-value references are moved and shared_ptr
 * are shared.
 *
 * \tparam T The target type
 * \tparam S The input type
 * \param s An l-value or r-value reference, or shared_ptr of type S
 * \returns Shared pointer of the base class type.
 */
template <class T, class S, typename = Impl::DisablePointer<S>>
auto wrap_own_share(S&& s) {
  return Impl::wrap_own_share<T>(std::forward<S>(s));
}

} // namespace Functions
} // namespace Dune
