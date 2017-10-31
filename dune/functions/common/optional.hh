// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_OPTIONAL_HH
#define DUNE_FUNCTIONS_COMMON_OPTIONAL_HH

#include <dune/common/std/optional.hh>

namespace Dune {
namespace Functions {
  template< class T >
  using Optional = Std::optional< T >;
}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_OPTIONAL_HH
