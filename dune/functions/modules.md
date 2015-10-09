<!-- -*- tab-width: 4; indent-tabs-mode: nil -*- -->
# Modules

@defgroup Functions
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

@defgroup FunctionInterface Polymorphic interface
@ingroup Functions
@brief Type-erasure based polymorphic interface

@defgroup FunctionImplementations Function implementations
@ingroup Functions
@brief Concrete function implementations







