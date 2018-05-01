// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ISTLVECTORBACKEND_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ISTLVECTORBACKEND_HH

#include <cstddef>
#include <utility>
#include <type_traits>

#include <dune/common/std/type_traits.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/functions/common/indexaccess.hh>
#include <dune/functions/functionspacebases/concepts.hh>


namespace Dune {
namespace Functions {



namespace Impl {



/*
 * \brief A wrapper providing multi-index access to vector entries
 *
 * The wrapped vector type should be an istl like random
 * access container providing operator[] and size() methods.
 * For classical containers this should support indices
 * of type std::size_t. For multi-type containers indices
 * of the form Dune::index_constant<i> should be supported
 * while size() should be a static constexpr method.
 *
 * When resolving multi-indices the backend appends indices
 * using operator[] as long as the result is not a scalar.
 * If this exhausts the digits of the multi-index, additional
 * zero`s are appended.
 *
 * \tparam V Type of the raw wrapper vector
 */
template<class V>
class ISTLVectorBackend
{

  // Template aliases for using detection idiom.
  template<class C>
  using dynamicIndexAccess_t = decltype(std::declval<C>()[0]);

  template<class C>
  using staticIndexAccess_t = decltype(std::declval<C>()[Dune::Indices::_0]);

  template<class C>
  using resizeMethod_t = decltype(std::declval<C>().resize(0));



  // Short cuts for feature detection
  template<class C>
  using hasDynamicIndexAccess = Dune::Std::is_detected<dynamicIndexAccess_t, std::remove_reference_t<C>>;

  template<class C>
  using hasStaticIndexAccess = Dune::Std::is_detected<staticIndexAccess_t, std::remove_reference_t<C>>;

  template<class C>
  using hasResizeMethod = Dune::Std::is_detected<resizeMethod_t, std::remove_reference_t<C>>;

  template<class C>
  using isDynamicVector = Dune::Std::is_detected<dynamicIndexAccess_t, std::remove_reference_t<C>>;

  template<class C>
  using isStaticVector = Dune::Std::bool_constant<
    Dune::Std::is_detected_v<staticIndexAccess_t, std::remove_reference_t<C>>
    and not Dune::Std::is_detected_v<dynamicIndexAccess_t, std::remove_reference_t<C>>>;

  template<class C>
  using isScalar = Dune::Std::bool_constant<not Dune::Std::is_detected_v<staticIndexAccess_t, std::remove_reference_t<C>>>;

  template<class C>
  using isVector = Dune::Std::bool_constant<Dune::Std::is_detected_v<staticIndexAccess_t, std::remove_reference_t<C>>>;



  // The method resolveMultiIndex() returns the entry of c
  // obtained by repeatedly calling [i_k] with all digits i_k
  // of the given multi-index i=(i_0, ..., i_m) starting from
  // k0=nextPosition. More precisely, it returns
  //
  //   c[i_k0][i_(k0+1)][i_(k0+2)]...[i_m][0]...[0]
  //
  // If the entry obtained after resolving all multi-index
  // digits is not scalar, i.e., if it still provides operator[]
  // access, then additional [0]'s are appended. This is needed
  // because Dune has no scalar type and Dune::FieldVector<K,1>
  // is commonly used for both, single component vectors and scalars.

  // This overload is used for dynamic vectors, i.e., for types
  // providing dynamic operator[] access.
  template<class C, class MultiIndex,
    std::enable_if_t<isDynamicVector<C>::value, int> = 0>
  static constexpr decltype(auto) resolveMultiIndex(C&& c, const MultiIndex& multiIndex, std::size_t nextPosition = 0)
  {
    std::size_t i = (nextPosition < multiIndex.size()) ? multiIndex[nextPosition] : 0;
    return resolveMultiIndex(c[i], multiIndex, nextPosition+1);
  }

  // This overload is used for static vectors, i.e., for types
  // providing (only) static operator[] access.
  template<class C, class MultiIndex,
    std::enable_if_t<isStaticVector<C>::value, int> = 0>
  static constexpr decltype(auto) resolveMultiIndex(C&& c, const MultiIndex& multiIndex, std::size_t nextPosition = 0)
  {
    std::size_t i = (nextPosition < multiIndex.size()) ? multiIndex[nextPosition] : 0;
    return hybridIndexAccess(c, i, [&] (auto&& ci) -> decltype(auto) {
      // Here we'd simply like to call resolveMultiIndex(...)
      // but gcc-5 does not find the name unless we explicitly
      // qualify the call although it should compile because
      // it's visible in our scope.
      return ISTLVectorBackend<V>::resolveMultiIndex(std::forward<decltype(ci)>(ci), multiIndex, nextPosition+1);
    });
  }

  // This overload is used for scalars, i.e., for types
  // that do not provide operator[] access.
  template<class C, class MultiIndex,
    std::enable_if_t<isScalar<C>::value, int> = 0>
  static constexpr decltype(auto) resolveMultiIndex(C&& c, const MultiIndex& multiIndex, std::size_t nextPosition = 0)
  {
    assert(nextPosition >= multiIndex.size());
    return std::forward<C>(c);
  }

  template<class... Args>
  static void forwardToResize(Args&&... args)
  {
    resize(std::forward<Args>(args)...);
  }


  template<class C, class SizeProvider,
    std::enable_if_t<hasResizeMethod<C>::value, int> = 0>
  static void resize(C&& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size==0)
    {
      // If size==0 this prefix refers to a single coefficient c.
      // But being in this overload means that c is not a scalar
      // because is has a resize() method. Since operator[] deliberately
      // supports implicit padding of multi-indices by as many
      // [0]'s as needed to resolve a scalar entry, we should also
      // except a non-scalar c here. However, this requires that
      // we silently believe that whatever size c already has is
      // intended by the user. The only exception is c.size()==0
      // which is not acceptable but we also cannot know the desired size.
      if (c.size()==0)
        DUNE_THROW(RangeError, "The vector entry v[" << prefix << "] should refer to a "
                                << "scalar coefficient, but is a dynamically sized vector of size==0");
      else
        // Accept non-zero sized coefficients to avoid that resize(basis)
        // fails for a vector that works with operator[] and already
        // has the appropriate size.
        return;
    }
    c.resize(size);
    prefix.push_back(0);
    for(std::size_t i=0; i<size; ++i)
    {
      prefix.back() = i;
      resize(c[i], sizeProvider, prefix);
    }
  }

  template<class C, class SizeProvider,
    std::enable_if_t<not hasResizeMethod<C>::value, int> = 0,
    std::enable_if_t<isVector<C>::value, int> = 0>
  static void resize(C&& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    // If size == 0 there's nothing to do:
    // We can't resize c and it's already
    // large enough anyway.
    if (size == 0)
      return;

    // If size>0 but c does not have the appropriate
    // size we throw an exception.
    //
    // We could perhaps also allow c.size()>size.
    // But then looping the loop below get's complicated:
    // We're not allowed to loop until c.size(). But
    // we also cannot use size for termination,
    // because this fails if c is a static vector.
    if (c.size() != size)
      DUNE_THROW(RangeError, "Can't resize non-resizable entry v[" << prefix << "] of size " << c.size() << " to size(" << prefix << ")=" << size);

    // Recursively resize all entries of c now.
    using namespace Dune::Hybrid;
    prefix.push_back(0);
    forEach(integralRange(Hybrid::size(c)), [&](auto&& i) {
      prefix.back() = i;
      // Here we'd simply like to call resize(c[i], sizeProvider, prefix);
      // but even gcc-7 does not except this bus reports
      // "error: ‘this’ was not captured for this lambda function"
      // although there's no 'this' because we're in a static method.
      // Bypassing this by and additional method that does perfect
      // forwarding allows to workaround this.
      ISTLVectorBackend<V>::forwardToResize(c[i], sizeProvider, prefix);
    });
  }

  template<class C, class SizeProvider,
    std::enable_if_t<not hasResizeMethod<C>::value, int> = 0,
    std::enable_if_t<isScalar<C>::value, int> = 0>
  static void resize(C&& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size != 0)
      DUNE_THROW(RangeError, "Can't resize scalar vector entry v[" << prefix << "] to size(" << prefix << ")=" << size);
  }

public:

  using Vector = V;

  ISTLVectorBackend(Vector& vector) :
    vector_(&vector)
  {}

  template<class SizeProvider>
  void resize(const SizeProvider& sizeProvider)
  {
    auto prefix = typename SizeProvider::SizePrefix();
    prefix.resize(0);
    resize(*vector_, sizeProvider, prefix);
  }

  template<class MultiIndex>
  decltype(auto) operator[](const MultiIndex& index) const
  {
    return resolveMultiIndex(*vector_, index);
  }

  template<class MultiIndex>
  decltype(auto) operator[](const MultiIndex& index)
  {
    return resolveMultiIndex(*vector_, index);
  }

  template<typename T>
  void operator= (const T& other)
  {
    vector() = other;
  }

  template<typename T>
  void operator= (const ISTLVectorBackend<T>& other)
  {
    vector() = other.vector();
  }

  const Vector& vector() const
  {
    return *vector_;
  }

  Vector& vector()
  {
    return *vector_;
  }

private:

  Vector* vector_;
};

} // end namespace Impl



/**
 * \brief Return a vector backend wrapping non-const ISTL like containers
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * The returned object implements the VectorBackend concept and
 * can be used for all dune-functions
 * utilities requiring a coefficient vector (e.g. interpolate()
 * and DiscreteGlobalBasisFunction). It essentially provides
 * operator[] access using multi-indices and a recursive
 * resize(GlobalBasis) method for adjusting the size to a
 * given GlobalBasis.
 *
 * Additionally to the VectorBackend interface, provides access
 * to the wrapped vector using the method vector() and forwards
 * all assignments to the underlying wrapped vector.
 *
 * The wrapped vector type should be a nested ISTL like random
 * access container providing operator[] and size() methods.
 * For classical containers this should support indices
 * of type std::size_t. For multi-type containers indices
 * of the form Dune::index_constant<i> should be supported
 * while size() should be a static constexpr method.
 *
 * When accessing the vector with a multi-index the backend
 * appends multi-index digits using operator[] as long as the
 * result is not a scalar. If this exhausts all digits of the
 * multi-index, additional zero`s are appended.
 *
 * \tparam V Type of the raw wrapper vector
 */
template<class Vector>
auto istlVectorBackend(Vector& v)
{
  return Impl::ISTLVectorBackend<Vector>(v);
}



/**
 * \brief Return a vector backend wrapping const ISTL like containers
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * The returned object implements the VectorBackend concept and
 * can be used for all dune-functions
 * utilities requiring a coefficient vector (e.g. interpolate()
 * and DiscreteGlobalBasisFunction. It essentially provides
 * operator[] access using multi-indices and a recursive
 * resize(GlobalBasis) method for adjusting the size to a given GlobalBasis.
 *
 * Additionally to the VectorBackend interface, provides access
 * to the wrapped vector using the method vector().
 *
 * The wrapped vector type should be a nested ISTL like random
 * access container providing operator[] and size() methods.
 * For classical containers this should support indices
 * of type std::size_t. For multi-type containers indices
 * of the form Dune::index_constant<i> should be supported
 * while size() should be a static constexpr method.
 *
 * When accessing the vector with a multi-index the backend
 * appends multi-index digits using operator[] as long as the
 * result is not a scalar. If this exhausts all digits of the
 * multi-index, additional zero`s are appended.
 *
 * \tparam V Type of the raw wrapper vector
 */
template<class Vector>
auto istlVectorBackend(const Vector& v)
{
  return Impl::ISTLVectorBackend<const Vector>(v);
}



} // namespace Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ISTLVECTORBACKEND_HH
