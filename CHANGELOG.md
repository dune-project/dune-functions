# Dune-functions changes

Any version of dune-functions is supposed to be compatible with the
correponding version of the Dune core modules.

## Master (will become release 2.7)

- `FlatVectorBackend` is now officially an implementation detail and thus moved
  to the namespace `Impl::`. The header `flatvectorbackend.hh` was removed.
  As a replacement the new free function `flatVectorView(c)` create a view
  object providing `operator[]` and `size()` methods for flat-vector-like
  access to the underlying container `c`.
- The `BasisBuilder` namespace has been renamed to `BasisFactory`.
  The old name still exists and should work as before, but its use
  is discouraged.

## Release 2.6

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
- An implementation of a global Brezzi-Douglas-Marini-Basis was added.  Thanks to Jakub Both
  for the code.
- `Dune::Functions::TupleVector` was deprecated. Use `Dune::TupleVector` from dune-common instead.

## Release 2.5

TODO...



