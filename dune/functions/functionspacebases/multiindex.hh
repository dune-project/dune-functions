// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MULTIINDEX_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MULTIINDEX_HH

#include <array>
#include <type_traits>

#include <dune/common/reservedvector.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {
namespace Impl {

template <class PreBasis>
struct MultiIndexSize
{
  static const std::size_t min = 1;
  static const std::size_t max = 1;
};

template <class PreBasis>
struct MultiIndexType
{
  static const std::size_t min = MultiIndexSize<PreBasis>::min;
  static const std::size_t max = MultiIndexSize<PreBasis>::max;
  using type = std::conditional_t<(max == 1), FlatMultiIndex<std::size_t>,
               std::conditional_t<(min == max), Dune::ReservedVector<std::size_t, max>, // std::array<std::size_t, max>,
                                                Dune::ReservedVector<std::size_t, max>>>;
};

} // end namespace Impl


template <class PreBasis>
using MultiIndexType_t = typename Impl::MultiIndexType<PreBasis>::type;

template <class PreBasis>
using SizePrefixType_t = Dune::ReservedVector<std::size_t, Impl::MultiIndexSize<PreBasis>::max>;


template <class MultiIndex>
void multiIndexPushFront(MultiIndex& M, std::size_t M0 = 0)
{
  M.resize(M.size()+1);
  for(std::size_t i=M.size()-1; i>0; --i)
    M[i] = M[i-1];
  M[0] = M0;
}

template <std::size_t N>
void multiIndexPushFront(std::array<std::size_t,N>& M, std::size_t M0 = 0)
{
  for(std::size_t i=M.size()-1; i>0; --i)
    M[i] = M[i-1];
  M[0] = M0;
}

template <class MultiIndex>
void multiIndexPushBack(MultiIndex& M, std::size_t M0 = 0)
{
  M.push_back(M0);
}

template <std::size_t N>
void multiIndexPushBack(std::array<std::size_t,N>& M, std::size_t M0 = 0)
{
  assert(false);
  // do nothing
}

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MULTIINDEX_HH
