// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORBACKEND_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORBACKEND_HH


#include <dune/common/concept.hh>

#include <dune/functions/functionspacebases/concepts.hh>




namespace Dune {
namespace Functions {



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


} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATVECTORBACKEND_HH
