// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_CONCEPT_HH
#define DUNE_FUNCTIONS_COMMON_CONCEPT_HH

#include <type_traits>
#include <utility>

#include <dune/functions/common/typelist.hh>

namespace Dune {
namespace Functions {
namespace Concept {

namespace Imp // All types and functions in this namespace are implementation details
{

  // Base class to mark refined concepts
  // Please derive from the derived template Concept::Refines
  struct RefinedConcept {};

  template<class C>
  static constexpr bool isRefinedConcept()
  { return std::is_base_of<RefinedConcept, C>(); }


  template<class C, class... T>
  static constexpr bool modelsImp();

  // Here is the implementation of the concept checking.
  // The first two overloads do the magic for checking
  // if the requirements of a concept are satisfied.
  // The rest is just for checking base concepts in case
  // of refinement.

  // This overload is present if type substitution for
  // C::require(T...) is successful, i.e., if the T...
  // matches the requirement of C. In this case this
  // overload is selected because C* is a better match than void*.
  template<class C, class... T,
    decltype(std::declval<C>().require(std::declval<T>()...), 0) =0>
  static constexpr auto matchesRequirement(C*) -> std::true_type
  { return std::true_type(); }

  // If the above overload is ruled out by SFINAE because
  // the T... does snot match the requirements of C, then
  // this default overload drops in.
  template<class C, class... T>
  static constexpr auto matchesRequirement(void*) -> std::false_type
  { return std::false_type(); }

  // Wrap above check into nice constexpr function
  template<class C, class...T>
  static constexpr bool matchesRequirement()
  {
    return decltype(matchesRequirement<C, T...>(std::declval<C*>()))::value;
  }



  // An empty list C of concepts is always matched by T...
  template<class C, class...T,
    typename std::enable_if< isEmptyTypeList<C>(), int>::type =0>
  static constexpr bool modelsConceptList()
  {
    return true;
  }

  // A nonempty list C of concepts is modeled
  // by T...  if it models the concept C::Head
  // and Concepts in C::Tail.
  template<class C, class...T,
    typename std::enable_if< not(isEmptyTypeList<C>()), int>::type =0>
  static constexpr bool modelsConceptList()
  {
    return modelsImp<typename C::Head, T...>() and modelsConceptList<typename C::Tail, T...>();
  }



  // If C is an unrefined concept, then T... models C
  // if it matches the requirement of C.
  template<class C, class... T,
    typename std::enable_if< not(isTypeList<C>()) and not(isRefinedConcept<C>()), int>::type=0>
  static constexpr bool modelsConcept()
  {
    return matchesRequirement<C, T...>();
  }

  // If C is a refined concept, then T... models C
  // if it matches the requirement of C and of
  // all base concepts.
  template<class C, class... T,
    typename std::enable_if< not(isTypeList<C>()) and isRefinedConcept<C>(), int>::type=0>
  static constexpr bool modelsConcept()
  {
    return matchesRequirement<C, T...>() and modelsConceptList<typename C::BaseConceptList, T...>();
  }

  // If C is a list of concepts, then T... models C
  // if matches the requirement of all concepts
  // in the list.
  template<class C, class... T,
    typename std::enable_if< isTypeList<C>(), int>::type=0>
  static constexpr bool modelsConcept()
  {
    return modelsConceptList<C, T...>();
  }

  // Check if T... models the concept or TypeList C
  template<class C, class... T>
  static constexpr bool modelsImp()
  {
    return modelsConcept<C, T...>();
  }
} // namespace Imp



// Check if T... models the concept or TypeList C
template<class C, class... T>
static constexpr bool models()
{
  return Imp::modelsImp<C, T...>();
}

// If you want to require
// A and B in a new concept C you
// should derive it from Refines<A,B>.
// By introducing the base class RefinedConcept
// we tag all refined concepts for recognition
// in overloads.
template<class... Base>
struct Refines : Imp::RefinedConcept
{
  typedef TypeList<Base...> BaseConceptList;
};



// Helper function for use in concept definitions.
// If the passed value b is not true, the concept will to be satisfied.
template<bool b, typename std::enable_if<b, int>::type = 0>
static constexpr bool requireTrue()
{
  return true;
}



}}} // namespace Dune::Functions::Concept




#endif // DUNE_FUNCTIONS_COMMON_CONCEPT_HH
