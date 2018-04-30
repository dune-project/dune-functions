// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORVIEW_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORVIEW_HH


#include <dune/common/concept.hh>

#include <dune/functions/functionspacebases/concepts.hh>




namespace Dune {
namespace Functions {
namespace Impl {



template<class V>
struct FlatVectorBackend
{

  template<class VV, class Index,
    typename std::enable_if< models<Concept::HasIndexAccess, VV, Index>(), int>::type = 0>
  static auto getEntry(VV&& v, const Index& i)
    ->decltype(v[i])
  {
    return v[i];
  }

  template<class VV, class Index,
    typename std::enable_if< not models<Concept::HasIndexAccess, VV, Index>(), int>::type = 0>
  static auto getEntry(VV&& v, const Index& i)
    ->decltype(v)
  {
    return std::forward<VV>(v);
  }

  template<class VV,
    typename std::enable_if< models<Concept::HasSizeMethod, VV>(), int>::type = 0>
  static auto size(VV&& v)
    ->decltype(v.size())
  {
    return v.size();
  }

  template<class VV,
    typename std::enable_if< not models<Concept::HasSizeMethod, VV>(), int>::type = 0>
  static std::size_t size(VV&& v)
  {
    return 1;
  }

};





template<class K, int n, int m>
struct FlatVectorBackend<typename Dune::FieldMatrix<K, n, m> >
{

  template<class VV, class Index>
  static auto getEntry(VV&& v, const Index& i) -> decltype(v[i/m][i%m])
  {
    return v[i/m][i%m];
  }

  template<class VV>
  static int size(VV&& v)
  {
    return n*m;
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
