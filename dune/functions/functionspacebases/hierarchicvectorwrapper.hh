// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORWRAPPER_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORWRAPPER_HH

#ifndef DUNE_FUNCTIONS_HIERARCHICVECTORWRAPPER_TEST_NO_DEPRECATION
#warning The header dune/functions/functionspacebases/hierarchicvectorwrapper.hh is deprecated and will be removed after release 2.10.
#endif

#include <dune/common/concept.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/common/indexaccess.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/concepts.hh>


namespace Dune {
namespace Functions {



namespace Imp {

  // Construct default coefficient type from vector and multiindex type
  // This requires that MultiIndex has a static size. Otherwise the
  // vector type itself is returned.
  template<class V, class MultiIndex>
  struct CoefficientType
  {
    template<class E, std::size_t size>
    struct DefaultCoefficientTypeHelper
    {
      using E0 = decltype(std::declval<E>()[Dune::Indices::_0]);
      using type = typename DefaultCoefficientTypeHelper<E0, size-1>::type;
    };

    template<class E>
    struct DefaultCoefficientTypeHelper<E, 0>
    {
      using type = E;
    };

    using type = typename DefaultCoefficientTypeHelper<V, StaticSizeOrZero<MultiIndex>::value>::type;
  };



  // This tag class is used as Coefficient template parameter
  // for HierarchicVectorWrapper if the coefficient type should
  // be deduced.
  struct DeducedCoefficientTag {};

} // namespace Imp



/**
 * \brief A wrapper providing multiindex access to vector entries
 *
 * The coefficient type should be a type such that the coefficients
 * entries for each global basis function can be cast to this type.
 * This is necessary because the wrapper cannot determine this type
 * automatically for multi-type containers and non-uniform indices.
 * The reason for this is, that the multi-index type will then be
 * dynamically sized such that the index depth cannot statically
 * be determined from the multi-indices. However, the compiler needs
 * a fixed termination criterion for instantiation of recursive
 * functions.
 *
 * If no coefficient type is given, the wrapper tries to determine
 * the coefficient type on its own assuming that the multi-indices
 * have fixed size.
 *
 * \deprecated This class is deprecated. Use istlVectorBackend(V) instead.
 *
 * \tparam V Type of the raw wrapper vector
 * \tparam CO Coefficient type
 */
template<class V, class CO=Imp::DeducedCoefficientTag>
class
[[deprecated("HierarchicVectorWrapper is deprecated and will be removed after release 2.10.")]]
HierarchicVectorWrapper
{
  template<class MultiIndex>
  using Coefficient = std::conditional_t< std::is_same_v<Imp::DeducedCoefficientTag,CO> and HasStaticSize_v<MultiIndex>,
            typename Imp::CoefficientType<V, MultiIndex>::type,
            CO
            >;


  using size_type = std::size_t;

  template<class C, class SizeProvider,
    std::enable_if_t< not models<Concept::HasResize, C>(), int> = 0,
    std::enable_if_t< not models<Concept::HasSizeMethod, C>(), int> = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size != 0)
      DUNE_THROW(RangeError, "Can't resize scalar vector entry v[" << prefix << "] to size(" << prefix << ")=" << size);
  }

  struct StaticResizeHelper
  {
    template<class I, class C, class SizeProvider>
    static void apply(I&& i, C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
    {
      prefix.back() = i;
      resizeHelper(c[i], sizeProvider, prefix);
    }
  };

  template<class C, class SizeProvider,
    std::enable_if_t< not models<Concept::HasResize, C>(), int> = 0,
    std::enable_if_t< models<Concept::HasSizeMethod, C>(), int> = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size == 0)
      return;

    if (c.size() != size)
      DUNE_THROW(RangeError, "Can't resize statically sized vector entry v[" << prefix << "] of size " << c.size() << " to size(" << prefix << ")=" << size);

    using namespace Dune::Hybrid;
    prefix.push_back(0);
    forEach(integralRange(Hybrid::size(c)), [&](auto&& i) {
        StaticResizeHelper::apply(i, c, sizeProvider, prefix);
      });
  }

  template<class C, class SizeProvider,
    std::enable_if_t< models<Concept::HasResize, C>(), int> = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size==0)
    {
      if (c.size()==0)
        DUNE_THROW(RangeError, "Can't resize dynamically sized vector entry v[" << prefix << "]. Its size is 0 but the target size is unknown due to size(" << prefix << ")=0.");
      else
        return;
    }

    c.resize(size);
    prefix.push_back(0);
    for(std::size_t i=0; i<size; ++i)
    {
      prefix.back() = i;
      resizeHelper(c[i], sizeProvider, prefix);
    }
  }



public:

  using Vector = V;

  template<class MultiIndex>
  using Entry = Coefficient<MultiIndex>;

  HierarchicVectorWrapper(Vector& vector) :
    vector_(&vector)
  {}

  template<class SizeProvider>
  void resize(const SizeProvider& sizeProvider)
  {
    typename SizeProvider::SizePrefix prefix;
    prefix.resize(0);
    resizeHelper(*vector_, sizeProvider, prefix);
  }

  template<class MultiIndex>
  const Entry<MultiIndex>& operator[](const MultiIndex& index) const
  {
      static_assert(not std::is_same_v<Imp::DeducedCoefficientTag,Entry<MultiIndex>>, "Coefficient type for HierarchicVectorWrapper and given multi-index type cannot be determined automatically!");
      return hybridMultiIndexAccess<const Entry<MultiIndex>&>(*vector_, index);
  }

  template<class MultiIndex>
  Entry<MultiIndex>& operator[](const MultiIndex& index)
  {
      static_assert(not std::is_same_v<Imp::DeducedCoefficientTag,Entry<MultiIndex>>, "Coefficient type for HierarchicVectorWrapper and given multi-index type cannot be determined automatically!");
      return hybridMultiIndexAccess<Entry<MultiIndex>&>(*vector_, index);
  }

  template<class MultiIndex>
  const Entry<MultiIndex>& operator()(const MultiIndex& index) const
  {
      static_assert(not std::is_same_v<Imp::DeducedCoefficientTag,Entry<MultiIndex>>, "Coefficient type for HierarchicVectorWrapper and given multi-index type cannot be determined automatically!");
      return (*this)[index];
  }

  template<class MultiIndex>
  Entry<MultiIndex>& operator()(const MultiIndex& index)
  {
      static_assert(not std::is_same_v<Imp::DeducedCoefficientTag,Entry<MultiIndex>>, "Coefficient type for HierarchicVectorWrapper and given multi-index type cannot be determined automatically!");
      return (*this)[index];
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



/**
 * \deprecated This function is deprecated. Use istlVectorBackend() instead.
 */
template<class V>
HierarchicVectorWrapper< V > hierarchicVector(V& v)
{
  return HierarchicVectorWrapper<V>(v);
}



/**
 * \deprecated This function is deprecated. Use istlVectorBackend() instead.
 */
template<class MultiIndex, class V,
    std::enable_if_t< models<Concept::HasIndexAccess, V, MultiIndex>(), int> = 0>
V& makeHierarchicVectorForMultiIndex(V& v)
{
  return v;
}



/**
 * \deprecated This function is deprecated. Use istlVectorBackend() instead.
 */
template<class MultiIndex, class V,
    std::enable_if_t< not models<Concept::HasIndexAccess, V, MultiIndex>(), int> = 0>
HierarchicVectorWrapper< V > makeHierarchicVectorForMultiIndex(V& v)
{
  return HierarchicVectorWrapper<V>(v);
}



} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORWRAPPER_HH
