// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORBACKEND_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORBACKEND_HH


#include <dune/functions/common/concept.hh>




namespace Dune {
namespace Functions {



namespace Concept {

struct HasIndexAcces
{
  template<class C, class I>
  auto require(C&& c, I&& i) -> decltype(
    c[i]
  );
};

struct HasSizeMethod
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.size()
  );
};

} // namespace Dune::Functions::Concept




template<class V>
struct FlatVectorBackend
{

  template<class VV, class Index,
    typename std::enable_if< Concept::models<Concept::HasIndexAcces, VV, Index>(), int>::type = 0>
  static auto getEntry(VV&& v, const Index& i)
    ->decltype(v[i])
  {
    return v[i];
  }

  template<class VV, class Index,
    typename std::enable_if< not Concept::models<Concept::HasIndexAcces, VV, Index>(), int>::type = 0>
  static auto getEntry(VV&& v, const Index& i)
    ->decltype(v)
  {
    return std::forward<VV>(v);
  }

  template<class VV,
    typename std::enable_if< Concept::models<Concept::HasSizeMethod, VV>(), int>::type = 0>
  static auto size(VV&& v)
    ->decltype(v.size())
  {
    return v.size();
  }

  template<class VV,
    typename std::enable_if< not Concept::models<Concept::HasSizeMethod, VV>(), int>::type = 0>
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


} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORBACKEND_HH
