// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORWRAPPER_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORWRAPPER_HH


#include <dune/functions/common/concept.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/staticforloop.hh>
#include <dune/functions/functionspacebases/concepts.hh>


namespace Dune {
namespace Functions {



/**
 * \brief A wrapper providing multiindex acces to vector entries
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
 * \tparam V Type of the raw wrapper vector
 * \tparam CO Coefficient type
 */
template<class V, class CO=void>
class HierarchicVectorWrapper
{
  using Coefficient = CO;
  using size_type = std::size_t;

  template<class C, class SizeProvider,
    typename std::enable_if< not Concept::models<Concept::HasResize, C>(), int>::type = 0,
    typename std::enable_if< not Concept::models<Concept::HasSizeMethod, C>(), int>::type = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size != 0)
      DUNE_THROW(RangeError, "Can't resize scalar vector entry v[" << prefix << "] to size(" << prefix << ")=" << size);
  }

  struct StaticResizeHelper
  {
    template<class I, class C, class SizeProvider>
    void operator()(I&& i, C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
    {
      prefix.back() = i;
      resizeHelper(c[i], sizeProvider, prefix);
    }
  };

  template<class C, class SizeProvider,
    typename std::enable_if< not Concept::models<Concept::HasResize, C>(), int>::type = 0,
    typename std::enable_if< Concept::models<Concept::HasSizeMethod, C>(), int>::type = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size == 0)
      return;

    if (c.size() != size)
      DUNE_THROW(RangeError, "Can't resize statically sized vector entry v[" << prefix << "] of size " << c.size() << " to size(" << prefix << ")=" << size);

    prefix.push_back(0);
    staticForLoop<0, StaticSize<C>::value>(StaticResizeHelper(), c, sizeProvider, prefix);
  }

  template<class C, class SizeProvider,
    typename std::enable_if< Concept::models<Concept::HasResize, C>(), int>::type = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size==0)
    {
      if (c.size()==0)
        DUNE_THROW(RangeError, "Can't resize dynamically sized vector entry v[" << prefix << "]. It's size is 0 but the target size is unknown due to size(" << prefix << ")=0.");
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




  template<size_type i, class E, class MultiIndex, size_type ti,
      typename std::enable_if< (ti < StaticSize<E>::value), int>::type = 0>
  static auto dynamicToStaticIndexAccess(E&& e, const MultiIndex& multiIndex, const Dune::TypeTree::index_constant<ti>& tryIndex)
    -> Coefficient&
  {
    if (multiIndex[i] == tryIndex)
      return getEntry<i+1>(e[tryIndex], multiIndex);
    return dynamicToStaticIndexAccess<i>(std::forward<E>(e), multiIndex, Dune::TypeTree::index_constant<ti+1>());
  }

  template<size_type i, class E, class MultiIndex>
  static auto dynamicToStaticIndexAccess(E&& e, const MultiIndex& multiIndex, const Dune::TypeTree::index_constant<StaticSize<E>::value>& tryIndex)
    -> Coefficient&
  {
    DUNE_THROW(RangeError, "Can't access e[" << StaticSize<E>::value << "] for statically sized container of size " << StaticSize<E>::value);
    return getEntry<i+1>(e[Dune::TypeTree::Indices::_0], multiIndex);
  }

  // The getEntry functions recursively apply operator[] to the given object e.
  // Each call of getEntry<i>(e, multiIndex) computes the multiIndex[i]-th entry
  // of e, calls getEntry<i+1>(e[multiIndex[i]], multiIndex) and returns the result.

  // This is the overload for dynamically sized e.
  template<size_type i, class E, class MultiIndex,
      typename std::enable_if< not std::is_convertible<typename std::decay<E>::type, Coefficient>::value, int>::type = 0,
      typename std::enable_if< not HasStaticSize<E>::value, int>::type = 0>
  static auto getEntry(E&& e, const MultiIndex& multiIndex)
    -> Coefficient&
  {
    return getEntry<i+1>(e[multiIndex[i]], multiIndex);
  }

  // This is the overload for statically sized e. Since e may be a multi-type container,
  // we can in general not use e[multiIndex[i]]. Instead we have to do a dynamic->static
  // mapping of multiIndex[i] which is achieved by dynamicToStaticIndexAccess.
  template<size_type i, class E, class MultiIndex,
      typename std::enable_if< not std::is_convertible<typename std::decay<E>::type, Coefficient>::value, int>::type = 0,
      typename std::enable_if< HasStaticSize<E>::value, int>::type = 0>
  static auto getEntry(E&& e, const MultiIndex& multiIndex)
    -> Coefficient&
  {
    using namespace Dune::TypeTree::Indices;
    return dynamicToStaticIndexAccess<i>(std::forward<E>(e), multiIndex, Dune::TypeTree::Indices::_0);
  }

  // This is the overload for the coefficient type. If e can be cast to the Coefficient type
  // we must assume that we have already addressed the coefficient and stop the index access
  // recursion. Notice that a check if i has reached the size of multiIndex cannot be used
  // as termination criterion in general, because multiIndex can be dynamically sized. This
  // is especially true, if the multi-indices are non-uniform.
  template<size_type i, class E, class MultiIndex,
      typename std::enable_if< std::is_convertible<typename std::decay<E>::type, Coefficient>::value, int>::type = 0>
  static auto getEntry(E&& e, const MultiIndex& multiIndex)
    -> Coefficient&
  {
    return std::forward<E>(e);
  }

public:

  using Vector = V;

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
  auto operator[](const MultiIndex& index) const
      ->decltype(getEntry<0>(std::declval<Vector>(), index))
  {
    return getEntry<0>(*vector_, index);
  }

  template<class MultiIndex>
  auto operator[](const MultiIndex& index)
      ->decltype(getEntry<0>(std::declval<Vector>(), index))
  {
    return getEntry<0>(*vector_, index);
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



template<class V>
class HierarchicVectorWrapper<V, void>
{

  template<class C, class SizeProvider,
    typename std::enable_if< not Concept::models<Concept::HasResize, C>(), int>::type = 0,
    typename std::enable_if< not Concept::models<Concept::HasSizeMethod, C>(), int>::type = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size != 0)
      DUNE_THROW(RangeError, "Can't resize scalar vector entry v[" << prefix << "] to size(" << prefix << ")=" << size);
  }

  struct StaticResizeHelper
  {
    template<class I, class C, class SizeProvider>
    void operator()(I&& i, C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
    {
      prefix.back() = i;
      resizeHelper(c[i], sizeProvider, prefix);
    }
  };

  template<class C, class SizeProvider,
    typename std::enable_if< not Concept::models<Concept::HasResize, C>(), int>::type = 0,
    typename std::enable_if< Concept::models<Concept::HasSizeMethod, C>(), int>::type = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size == 0)
      return;

    if (c.size() != size)
      DUNE_THROW(RangeError, "Can't resize statically sized vector entry v[" << prefix << "] of size " << c.size() << " to size(" << prefix << ")=" << size);

    prefix.push_back(0);
    staticForLoop<0, StaticSize<C>::value>(StaticResizeHelper(), c, sizeProvider, prefix);
  }

  template<class C, class SizeProvider,
    typename std::enable_if< Concept::models<Concept::HasResize, C>(), int>::type = 0>
  static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
  {
    auto size = sizeProvider.size(prefix);
    if (size==0)
    {
      if (c.size()==0)
        DUNE_THROW(RangeError, "Can't resize dynamically sized vector entry v[" << prefix << "]. It's size is 0 but the target size is unknown due to size(" << prefix << ")=0.");
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

  template<int start, int end>
  struct GetEntryHelper
  {
    template<class T, class MultiIndex,
      typename std::enable_if< Concept::models<Concept::HasResize, typename std::decay<T>::type>(), int>::type = 0>
    static auto getEntry(T&& t, MultiIndex&& index)
      -> decltype(GetEntryHelper<start+1,end>::getEntry(t[index[start]], index))
    {
      return GetEntryHelper<start+1,end>::getEntry(t[index[start]], index);
    }

    template<class T, class MultiIndex,
      typename std::enable_if< not Concept::models<Concept::HasResize, typename std::decay<T>::type>(), int>::type = 0>
    static auto getEntry(T&& t, MultiIndex&& index)
      -> decltype(GetEntryHelper<start+1,end>::getEntry(t[Dune::TypeTree::Indices::_0], index))
    {
      using namespace Dune::TypeTree::Indices;
      static const int size = StaticSize<typename std::decay<T>::type>::value;
      return getEntry(std::forward<T>(t), index, Dune::TypeTree::index_constant<size-1>());
    }

    template<class T, class MultiIndex, std::size_t i,
      typename std::enable_if< (i>0), int>::type = 0>
    static auto getEntry(T&& t, MultiIndex&& index, const Dune::TypeTree::index_constant<i>& static_i)
      -> decltype(GetEntryHelper<start+1,end>::getEntry(t[Dune::TypeTree::Indices::_0], index))
    {
      if (index[start]==i)
        return GetEntryHelper<start+1,end>::getEntry(t[static_i], index);
      return GetEntryHelper<start,end>::getEntry(std::forward<T>(t), index, Dune::TypeTree::index_constant<i-1>());
    }

    template<class T, class MultiIndex>
    static auto getEntry(T&& t, MultiIndex&& index, const Dune::TypeTree::index_constant<0>& static_i)
      -> decltype(GetEntryHelper<start+1,end>::getEntry(t[static_i], index))
    {
      return GetEntryHelper<start+1,end>::getEntry(t[static_i], index);
    }

  };

  template<int start>
  struct GetEntryHelper<start, start>
  {
    template<class T, class MultiIndex>
    static auto getEntry(T&& t, MultiIndex&& index)
      ->decltype(std::forward<T>(t))
    {
      return std::forward<T>(t);
    }
  };

public:

  using Vector = V;

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

  template<class MultiIndex,
           typename std::enable_if< HasStaticSize<MultiIndex>::value, int>::type = 0>
  auto operator[](MultiIndex&& index) const
      ->decltype(GetEntryHelper<0, StaticSize<MultiIndex>::value>::getEntry(std::declval<Vector>(), index))
  {
    return GetEntryHelper<0, StaticSize<MultiIndex>::value>::getEntry(*vector_, index);
  }

  template<class MultiIndex,
           typename std::enable_if< HasStaticSize<MultiIndex>::value, int>::type = 0>
  auto operator[](MultiIndex&& index)
      ->decltype(GetEntryHelper<0, StaticSize<MultiIndex>::value>::getEntry(std::declval<Vector>(), index))
  {
    return GetEntryHelper<0, StaticSize<MultiIndex>::value>::getEntry(*vector_, index);
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



template<class V>
HierarchicVectorWrapper< V > hierarchicVector(V& v)
{
  return HierarchicVectorWrapper<V>(v);
}



template<class MultiIndex, class V,
    typename std::enable_if< Concept::models<Concept::HasIndexAccess, V, MultiIndex>(), int>::type = 0>
V& makeHierarchicVectorForMultiIndex(V& v)
{
  return v;
}



template<class MultiIndex, class V,
    typename std::enable_if< not Concept::models<Concept::HasIndexAccess, V, MultiIndex>(), int>::type = 0>
HierarchicVectorWrapper< V > makeHierarchicVectorForMultiIndex(V& v)
{
  return HierarchicVectorWrapper<V>(v);
}



} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORWRAPPER_HH
