// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MULTIINDEX_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MULTIINDEX_HH

#include <type_traits>

#include <dune/common/reservedvector.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune { namespace Functions
{
  template <class PreBasis>
  struct RequiredMultiIndexSize
  {
    static constexpr std::size_t value = 1;
  };


  namespace Impl
  {
    template <class PreBasis>
    struct MultiIndexType
    {
      static constexpr std::size_t N = RequiredMultiIndexSize<PreBasis>::value;
      using type = std::conditional_t<(N == 1),
        FlatMultiIndex<std::size_t>,
        Dune::ReservedVector<std::size_t, N>>;
    };
  }

  template <class PreBasis>
  using MultiIndexType_t = typename Impl::MultiIndexType<PreBasis>::type;

  template <class PreBasis>
  using SizePrefixType_t = Dune::ReservedVector<std::size_t, RequiredMultiIndexSize<PreBasis>::value>;


}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MULTIINDEX_HH
