# Dune-functions changes

Any version of dune-functions is supposed to be compatible with the
corresponding version of the Dune core modules.

## Master (will become release 2.8)

## Release 2.7

- The `LagrangeBasis` is extended by a template parameter to set the range type of
  the underlying LocalBasis. This parameter is also added to the basis factory.
  One can write `lagrange<k, float>()` to create a lagrange basis with compile-time
  order `k` and range type `float`. The range type defaults to `double` if
  nothing is given. The run-time order lagrange functions can be created by
  `lagrange<float>(k)`, correspondingly.
- The `LagrangeBasis` class can now be used with a run-time polynomial order.
  For this, the template parameter setting the order now has a default value of -1.
  If this default value is used, the class accepts an integer constructor
  argument with a run-time order.  As the decision whether to use a compile-time
  or run-time order is taken at compile-time there should be no efficiency
  decrease for the compile-time-order case.
  (Note: The implementation currently uses the `LagrangeFiniteElement`
  implementation of `dune-localfunctions`, which violates strict-aliasing rules.
  This may lead to surprising crashes and miscalculations when compiling
  with optimization.)

## Release 2.6

- The functionality provided by `LocalIndexSet`, namely the `indices`
  method and the `MultiIndex` typedef, have been merged into
  `LocalView`.  The global multi-indices are now precomputed and cached
  when binding the `LocalView`. `LocalIndexSet` as well as the `localIndexSet()`
  method of global bases have been deprecated.
- When using `DiscreteGlobalBasisFunction` and `interpolate` the
  function range is accessed hierarchically. I.e. a leaf node
  with tree path i0,...,in is mapped to the range entry `r[i0]...[in]`.
- If the vector passed to `DiscreteGlobalBasisFunction` or `interpolate()`
  does not provide the `VectorBackend` interface, it is automatically
  wrapped using `istlVectorBackend(vector)`.
- The new method `istlVectorBackend(vector)` creates a `VectorBackend`
  suitable for being used as coefficient vector for `interpolate()`
  and `DiscreteGlobalBasisFunction`. The underlying `vector` can
  be a nested vector composed by stl and dune random access containers.
  Notice that the only scalar coefficients per basis function are supported.
- The algorithm `forEachBoundaryDOF()` was added in a new header `boundarydofs.hh`.
  It allows to iterate over all DOFs of a given global basis associated to sub-entities
  located at the boundary of the domain.
- The class `SubEntityDOFs<GridView>` was added along with some `subEnitityDOFs`
  helper functions for creation in the new header `subentitydofs.hh`. After
  binding an object of this class to a local view and a sub-entity it can be
  used a range of local indices of DOFs associated to sub-sub-entities of this sub-entity.
- `FlatVectorBackend` is now officially an implementation detail and thus moved
  to the namespace `Impl::`. The header `flatvectorbackend.hh` was removed.
  As a replacement the new free function `flatVectorView(c)` create a view
  object providing `operator[]` and `size()` methods for flat-vector-like
  access to the underlying container `c`.
- The `BasisBuilder` namespace has been renamed to `BasisFactory`.
  The old name still exists and should work as before, but its use
  is discouraged.
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



