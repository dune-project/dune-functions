// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH


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



struct HasSizeMethod
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.size()
  );
};



struct HasIndexAccess
{
  template<class C, class I>
  auto require(C&& c, I&& i) -> decltype(
    c[i]
  );
};



} // namespace Dune::Functions::Concept
} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
