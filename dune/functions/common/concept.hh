// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_CONCEPT_HH
#define DUNE_FUNCTIONS_COMMON_CONCEPT_HH

#include <type_traits>
#include <utility>

#include <dune/common/prioritytag.hh>
#include <dune/functions/common/typelist.hh>

namespace Dune {
namespace Functions {
namespace Concept {

// #############################################################################
// # All types and functions following here are implementation details
// # for the models() function below.
// #############################################################################
namespace Imp
{

  // Base class to mark refined concepts
  // Please derive from the derived template Concept::Refines
  struct RefinedConcept {};

  template<class C>
  constexpr bool isRefinedConcept()
  { return std::is_base_of<RefinedConcept, C>(); }


  template<class C, class... T>
  constexpr bool modelsImp();

  // Here is the implementation of the concept checking.
  // The first two overloads do the magic for checking
  // if the requirements of a concept are satisfied.
  // The rest is just for checking base concepts in case
  // of refinement.

  // This overload is present if type substitution for
  // C::require(T...) is successful, i.e., if the T...
  // matches the requirement of C. In this case this
  // overload is selected because PriorityTag<1>
  // is a better match for PrioriryTag<42> than
  // PriorityTag<0> in the default overload.
  template<class C, class... T,
    decltype(std::declval<C>().require(std::declval<T>()...), 0) =0>
  constexpr std::true_type matchesRequirement(PriorityTag<1>)
  { return {}; }

  // If the above overload is ruled out by SFINAE because
  // the T... does not match the requirements of C, then
  // this default overload drops in.
  template<class C, class... T>
  constexpr std::false_type matchesRequirement(PriorityTag<0>)
  { return {}; }

  // Wrap above check into nice constexpr function
  template<class C, class...T>
  constexpr auto matchesRequirement()
    -> decltype(matchesRequirement<C, T...>(PriorityTag<42>()))
  { return {}; }



  // An empty list C of concepts is always matched by T...
  template<class C, class...T,
    typename std::enable_if< isEmptyTypeList<C>(), int>::type =0>
  constexpr std::true_type modelsConceptList()
  { return {}; }

  // A nonempty list C of concepts is modeled
  // by T...  if it models the concept C::Head
  // and Concepts in C::Tail.
  template<class C, class...T,
    typename std::enable_if< not(isEmptyTypeList<C>()), int>::type =0>
  constexpr auto modelsConceptList()
    -> std::integral_constant<bool, modelsImp<typename C::Head, T...>() and modelsConceptList<typename C::Tail, T...>()>
  { return {}; }

  // If C is an unrefined concept, then T... models C
  // if it matches the requirement of C.
  template<class C, class... T,
    typename std::enable_if< not(isTypeList<C>()) and not(isRefinedConcept<C>()), int>::type=0>
  constexpr auto modelsConcept()
    -> decltype(matchesRequirement<C, T...>())
  { return {}; }

  // If C is a refined concept, then T... models C
  // if it matches the requirement of C and of
  // all base concepts.
  template<class C, class... T,
    typename std::enable_if< not(isTypeList<C>()) and isRefinedConcept<C>(), int>::type=0>
  constexpr auto modelsConcept()
    -> std::integral_constant<bool, matchesRequirement<C, T...>() and modelsConceptList<typename C::BaseConceptList, T...>()>
  { return {}; }

  // If C is a list of concepts, then T... models C
  // if matches the requirement of all concepts
  // in the list.
  template<class C, class... T,
    typename std::enable_if< isTypeList<C>(), int>::type=0>
  constexpr auto modelsConcept()
    -> decltype(modelsConceptList<C, T...>())
  { return {}; }

  // Check if T... models the concept or TypeList C
  // Here we cannot use true_type/false_type as return
  // type because we need a recursion for checking base
  // concepts. In order to make the recursion compile we
  // need a forward declaration of the corresponding function,
  // which is not possible if the result is encoded in the
  // return type and not only the return value.
  template<class C, class... T>
  constexpr bool modelsImp()
  {
    return modelsConcept<C, T...>();
  }
} // namespace Imp



// #############################################################################
// # The method models() does the actual check if a type models a concept
// # using the implementation details above.
// #############################################################################

// Check if T... models the concept or TypeList C
template<class C, class... T>
constexpr bool models()
{
  return Imp::modelsImp<C, T...>();
}



// #############################################################################
// # All functions following here are implementation details for the
// # for the tupleEntriesModel() function below.
// #############################################################################
namespace Imp
{

  template<class C, class First>
  constexpr auto allModel()
    -> std::integral_constant<bool, Concept::models<C, First>()>
  { return {}; }

  template<class C, class First, class... Other>
  constexpr auto allModel()
    -> std::integral_constant<bool, Concept::models<C, First>() and allModel<C, Other...>()>
  { return {}; }

  template<class C, class... T>
  constexpr auto tupleEntriesModel(const std::tuple<T...>&)
    -> decltype(allModel<C, T...>())
  { return {}; }

}



// #############################################################################
// # The method tupleEntriesModel() does the actual check if the types in a tuple
// # model a concept using the implementation details above.
// #############################################################################

template<class C, class Tuple>
constexpr auto tupleEntriesModel()
  -> decltype(Imp::tupleEntriesModel<C>(std::declval<Tuple>()))
{
  return {};
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



// #############################################################################
// # The following require*() functions are just helpers that allow to
// # propagate a failed check as substitution failure. This is usefull
// # inside of a concept definition.
// #############################################################################

// Helper function for use in concept definitions.
// If the passed value b is not true, the concept will to be satisfied.
template<bool b, typename std::enable_if<b, int>::type = 0>
constexpr bool requireTrue()
{
  return true;
}

// Helper function for use in concept definitions.
template<class C, class... T, typename std::enable_if<models<C, T...>(), int>::type = 0>
constexpr bool requireConcept()
{
  return true;
}

// Helper function for use in concept definitions.
// This allows to avoid using decltype
template<class C, class... T, typename std::enable_if<models<C, T...>(), int>::type = 0>
constexpr bool requireConcept(T&&... t)
{
  return true;
}

// Helper function for use in concept definitions.
// This checks if the concept given as first type is modelled by all types in the tuple passed as argument
template<class C, class Tuple, typename std::enable_if<tupleEntriesModel<C, Tuple>(), int>::type = 0>
constexpr bool requireConceptForTupleEntries()
{
  return true;
}

// Helper function for use in concept definitions.
// If the first passed type is not convertible to the second, the concept will not be satisfied.
template<class From, class To,
  typename std::enable_if< std::is_convertible<From, To>::value, int>::type = 0>
constexpr bool requireConvertible()
{
  return true;
}

// Helper function for use in concept definitions.
// If passed argument is not convertible to the first passed type, the concept will not be satisfied.
template<class To, class From,
  typename std::enable_if< std::is_convertible<From, To>::value, int>::type = 0>
constexpr bool requireConvertible(const From&)
{
  return true;
}

// Helper function for use in concept definitions.
// This will always evaluate to true. If just allow
// to turn a type into an expression. The failure happens
// already during substitution for the type argument.
template<typename T>
constexpr bool requireType()
{
  return true;
}

// Helper function for use in concept definitions.
// If first passed type is not a base class of second type, the concept will not be satisfied.
template<class Base, class Derived,
  typename std::enable_if< std::is_base_of<Base, Derived>::value, int>::type = 0>
constexpr bool requireBaseOf()
{
  return true;
}

// Helper function for use in concept definitions.
// If first passed type is not a base class of first arguments type, the concept will not be satisfied.
template<class Base, class Derived,
  typename std::enable_if< std::is_base_of<Base, Derived>::value, int>::type = 0>
constexpr bool requireBaseOf(const Derived&)
{
  return true;
}

// Helper function for use in concept definitions.
// If the passed types are not the same, the concept will not be satisfied.
template<class A, class B,
  typename std::enable_if< std::is_same<A, B>::value, int>::type = 0>
constexpr bool requireSameType()
{
  return true;
}



}}} // namespace Dune::Functions::Concept




#endif // DUNE_FUNCTIONS_COMMON_CONCEPT_HH
