// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TYPEUTILITIES_HH
#define DUNE_FUNCTIONS_COMMON_TYPEUTILITIES_HH

#include <memory>

#include <dune/common/std/type_traits.hh>

namespace Dune {
namespace Functions {

namespace Impl {

template <class T>
struct is_shared_ptr : std::false_type {};

template <class T>
struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {};

template <class T>
constexpr auto isSharedPtr() {
  return Std::bool_constant<is_shared_ptr<std::decay_t<T>>::value>();
}

} // end namespace Impl

template <class T, class R = void>
using enableIfSharedPtr = std::enable_if_t<Impl::isSharedPtr<T>(), R>;

template <class T, class R = void>
using disableIfSharedPtr = std::enable_if_t<not Impl::is_shared_ptr<T>::value, R>;

namespace Impl {

template <class T, class Enable = void> struct PointerDecay;

template <class T>
struct PointerDecay<T, disableIfSharedPtr<T>> {
  using type = std::decay_t<T>;
};

template <class T>
struct PointerDecay<T, enableIfSharedPtr<T>> {
  using RawSharedPointer = std::decay_t<T>;
  using WrappedType = typename RawSharedPointer::element_type;
  using type = std::decay_t<WrappedType>;
};

} // namespace Impl

// decays a single layer of shared_ptr wrapping (if there is such)
// and the cv qualifiers from the remaining type
template <class TT>
using pointer_decay_t = typename Impl::PointerDecay<TT>::type;

}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_TYPEUTILITIES_HH
