<!-- -*- tab-width: 4; indent-tabs-mode: nil -*- -->
# Modules

@defgroup Functions Functions
@brief Interfaces and implementations of functions

Since interfaces in dune-functions rely on duck-typing
there are no base classes to ensure interfaces. Instead
of this each type can be checked, if it satisfies a certain
concept as defined in the Dune::Functions::Concept namespace.
Additionally there is a polymorphic interface consisting
of type-erasure wrappers in the spirit of std::function.

@defgroup FunctionConcepts Function concepts
@ingroup Functions
@brief Concept definitions for function interfaces

@defgroup FunctionInterface Function interface wrappers
@ingroup Functions
@brief Type-erasure based polymorphic wrappers

@defgroup FunctionImplementations Function implementations
@ingroup Functions
@brief Concrete function implementations

@defgroup FunctionUtility Function related utilities
@ingroup Functions
@brief Helper classes and functions related to functions



@defgroup FunctionSpaceBases Function space bases
@brief Interfaces and implementation of global bases for grid function space

@defgroup FunctionSpaceBasesConcepts Function space basis concepts
@ingroup FunctionSpaceBases

@defgroup FunctionSpaceBasesImplementations Function space basis implementations
@ingroup FunctionSpaceBases
@brief Implementations of the GlobalBasis concept

@defgroup FunctionSpaceBasesUtilities Function space basis related utilities
@ingroup FunctionSpaceBases



@defgroup Utility Utility
@brief Various helper classes and functions

@defgroup TypeErasure Utilities for type-erasure
@brief Helper classes for implementing type-erased interfaces
@ingroup Utility






