// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_BACKEND_CONCEPTS_HH
#define DUNE_FUNCTIONS_BACKEND_CONCEPTS_HH


#include <utility>

#include <dune/common/concept.hh>

namespace Dune {
namespace Functions {
namespace Concept {

using namespace Dune::Concept;


// Concept for a ConstVectorBackend
template<class GlobalBasis>
struct ConstVectorBackend
{
  template<class V>
  auto require(const V& v) -> decltype(
    v[std::declval<typename GlobalBasis::MultiIndex>()]
  );
};

// Concept for a VectorBackend
template<class GlobalBasis>
struct VectorBackend : Refines<ConstVectorBackend<GlobalBasis>>
{
  template<class V>
  auto require(const V& v) -> decltype(
    const_cast<V&>(v).resize(std::declval<const GlobalBasis&>()),
    const_cast<V&>(v)[std::declval<typename GlobalBasis::MultiIndex>()] = v[std::declval<typename GlobalBasis::MultiIndex>()]
  );
};

} // namespace Dune::Functions::Concept
} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_BACKEND_CONCEPTS_HH
