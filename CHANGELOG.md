dune-functions
========

This is the development version of dune-functions leading up to version 2.6.

It is supposed to be compatible with the Dune core modules of the same version.

Changes
=======

dune-functions 2.6-dev
------------

- Added an implementation of a Rannacher-Turek basis
- Add a set of unit tests for bases
- Extend the documentation of these bases
- Remove `DiscreteScalarGlobalBasisFunction`. This was replaced by `DiscreteGlobalBasisFunction`.
- Introduce inteface method `NodeIndexSet::indices(it)` that must write all global
  indices for the bound tree into the container pointed to by the iterator `it`
- Remove `NodeIndexSet::index(i)` for computation of individual global indices.
- `DefaultLocalIndexSet::bind(...)` now computes and caches all global indices
  while `DefaultLocalIndexSet::index(...)` is just a cheap lookup in the cache.
- `Dune::Functions::Optional` was removed. Please use `Dune::Std::optional` instead.
- Dune-functions now requires at least GCC 5 and uses C++14. As a consequence
  Debian 8 is no longer supported and tested.
- An implementation of a global Raviart-Thomas basis provided by Jakub Both was merged.
  This is still incomplete and considered experimental.

dune-functions 2.5
------------

TODO...



