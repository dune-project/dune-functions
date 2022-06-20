// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORVIEW_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORVIEW_HH


#include <array>

#include <dune/common/concept.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/functionspacebases/concepts.hh>




namespace Dune {
namespace Functions {
namespace Impl {


template<class V>
struct FlatVectorBackend
{

  template<class VV, class Index,
    typename std::enable_if< models<Concept::HasIndexAccess, VV, Index>(), int>::type = 0>
  static decltype(auto) getEntry(VV&& v, const Index& i)
  {
    return v[i];
  }

  template<class VV, class Index,
    typename std::enable_if< not models<Concept::HasIndexAccess, VV, Index>(), int>::type = 0>
  static decltype(auto) getEntry(VV&& v, const Index&)
  {
    return std::forward<VV>(v);
  }

  template<class VV,
    typename std::enable_if< models<Concept::HasSizeMethod, VV>(), int>::type = 0>
  static auto size(VV&& v)
  {
    return Dune::Hybrid::size(v);
  }

  template<class VV,
    typename std::enable_if< not models<Concept::HasSizeMethod, VV>(), int>::type = 0>
  static auto size(VV&&)
  {
    return Dune::index_constant<1>{};
  }
};




template<class K, int n, int m>
struct FlatVectorBackend<typename Dune::FieldMatrix<K, n, m> >
{

  template<class VV, class Index>
  static decltype(auto) getEntry(VV&& v, const Index& i)
  {
    return v[i/m][i%m];
  }

  template<class VV>
  static auto size(VV&& v)
  {
    return Dune::index_constant<n*m>{};
  }
};



template<class K, std::size_t n>
struct FlatVectorBackend< std::array<K, n> >
{

  template<class VV, class Index>
  static decltype(auto) getEntry(VV&& v, const Index& i)
  {
    const auto innerSize = decltype(FlatVectorBackend<K>::size(v[0]))::value;
    return FlatVectorBackend<K>::getEntry(v[i/innerSize], i%innerSize);
  }

  template<class VV>
  static auto size(VV&& v)
  {
    const auto innerSize = decltype(FlatVectorBackend<K>::size(v[0]))::value;
    return Dune::index_constant<n*innerSize>{};
  }

};




template<class T>
class FlatVectorView
{
  using Backend = FlatVectorBackend<std::decay_t<T>>;
public:
  FlatVectorView(T& t) :
    t_(&t)
  {}

  auto size() const
  {
    return Backend::size(*t_);
  }

  template<class Index>
  decltype(auto) operator[](const Index& i) const
  {
    return Backend::getEntry(*t_, i);
  }

  template<class Index>
  decltype(auto) operator[](const Index& i)
  {
    return Backend::getEntry(*t_, i);
  }

private:
  T* t_;
};


template<class T>
class FlatVectorView<T&&>
{
  using Backend = FlatVectorBackend<std::decay_t<T>>;
public:
  FlatVectorView(T&& t) :
    t_(std::move(t))
  {}

  auto size() const
  {
    return Backend::size(t_);
  }

  template<class Index>
  decltype(auto) operator[](const Index& i) const
  {
    return Backend::getEntry(t_, i);
  }

  template<class Index>
  decltype(auto) operator[](const Index& i)
  {
    return Backend::getEntry(t_, i);
  }

private:
  T t_;
};

} // namespace Impl



/**
 * \brief Create flat vector view of passed mutable container
 *
 * When passed a nested container, the resulting value is
 * a flat-vector-like view object. It provides an operator[]
 * method to access all entries of the underlying nested
 * container using flat consecutive indices and a size()
 * method to compute the corresponding total size.
 *
 * This method will create a view object storing a pointer
 * to the passed mutable container.
 */
template<class T>
auto flatVectorView(T& t)
{
  return Impl::FlatVectorView<T>(t);
}

/**
 * \brief Create flat vector view of passed const container
 *
 * When passed a nested container, the resulting value is
 * a flat-vector-like view object. It provides an operator[]
 * method to access all entries of the underlying nested
 * container using flat consecutive indices and a size()
 * method to compute the corresponding total size.
 *
 * This method will create a view object storing a pointer
 * to the passed const container.
 */
template<class T>
auto flatVectorView(const T& t)
{
  return Impl::FlatVectorView<const T>(t);
}

/**
 * \brief Create flat vector view of passed container temporary
 *
 * When passed a nested container, the resulting value is
 * a flat-vector-like view object. It provides an operator[]
 * method to access all entries of the underlying nested
 * container using flat consecutive indices and a size()
 * method to compute the corresponding total size.
 *
 * This method will create a 'view' object storing the
 * provided temporary container by value.
 */
template<class T>
auto flatVectorView(T&& t)
{
  return Impl::FlatVectorView<T&&>(std::move(t));
}


} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORVIEW_HH
