// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORBACKEND_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORBACKEND_HH


#include <dune/functions/common/concept.hh>


namespace Dune {
namespace Functions {



namespace Concept {

struct HasResize
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.resize(0)
  );
};

} // namespace Dune::Functions::Concept



struct HierarchicVectorBackend
{

    template<class C, class SizeProvider,
        typename std::enable_if< not Concept::models<Concept::HasResize, C>(), int>::type = 0>
    static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
    {}

    template<class C, class SizeProvider,
        typename std::enable_if< Concept::models<Concept::HasResize, C>(), int>::type = 0>
    static void resizeHelper(C& c, const SizeProvider& sizeProvider, typename SizeProvider::SizePrefix prefix)
    {
        auto size = sizeProvider.size(prefix);
        c.resize(size);
        prefix.push_back(0);
        for(std::size_t i=0; i<size; ++i)
        {
            prefix.back() = i;
            resizeHelper(c[i], sizeProvider, prefix);
        }
    }

    template<class V, class SizeProvider>
    static void resize(V& v, const SizeProvider& sizeProvider)
    {
        typename SizeProvider::SizePrefix prefix;
        prefix.resize(0);
        resizeHelper(v, sizeProvider, prefix);
    }


    template<int start, int end>
    struct GetEntryHelper
    {
        template<class T, class MultiIndex>
        static auto getEntry(T&& t, MultiIndex&& index)
            -> decltype(GetEntryHelper<start+1,end>::getEntry(t[index[start]], index))
        {
            return GetEntryHelper<start+1,end>::getEntry(t[index[start]], index);
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


    template<class VV, class MultiIndex>
    static auto getEntry(VV&& v, MultiIndex&& index)
        ->decltype(GetEntryHelper<0, std::tuple_size<typename std::decay<MultiIndex>::type>::value>::getEntry(v, index))
    {
        return GetEntryHelper<0, std::tuple_size<typename std::decay<MultiIndex>::type>::value>::getEntry(v, index);
    }

};



} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICVECTORBACKEND_HH
