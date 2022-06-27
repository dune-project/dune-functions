# Dune-functions changes

Any version of dune-functions is supposed to be compatible with the
corresponding version of the Dune core modules.

## Master (will become release 2.9)

- The `MultiIndex` used by `DefaultGlobalBasis` is now a `StaticMultiIndex`
  if it's size is known at compile time. Otherwise its a `ReservedVector`.
- The template class `StaticMultiIndex` for a statically sized multiindex was added.
  This essentially adds `operator<<` to `std::array` for writing to a stream
  and a cast to the first entry for `size()==0`.
  `FlatMultiIndex<T>` is now an alias for `StaticMultiIndex<T,1>`.
- The template class `OverflowArray` was added. It mostly behaves like `Dune::ReservedVector`
  but derives from a statically sized array base class. This allows to have temporary
  dynamic size but cast to the result to the statically sized base class.
- The new `FaceNormalGridFunction` implements the unit outer normal
  of the closest face within an element as a grid function.
- SubspaceBases will no longer be nested. Instead,
  `SubspaceBasis(SubspaceBases(rootBasis,innerTP),outerTP)`
  is resolved to a `SubspaceBases(rootBasis, tp)` where the tree path `tp`
  is obainted by joining `innerTP` and `outerTP`.
  This is implemented for the helper function `subspaceBasis()`
  and the newly added class template argument deduction guides.
- The way multi index types are generated for global bases has been refactored.
  This changes nothing for the user but gives more freedom in the construction
  of global bases.
  But it is a breaking change for the `PreBasis` and `PreBasisFactory` interface.
  If you maintain your own prebasis implementions, they have to be adapted.
  The necessary modifications are documented in detail here: https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/326.
- Support for `PreBasis` implementations not providing `PreBasis::indices()`
  has been removed.
- The added `ComposedGridFunction` implements the composition of a
  lambda function with one or several grid functions while itself
  providing the grid function interface.
- The new utility `resolveRef(t)` allows to seamlessly work with
  `std::reference_wrapper` by returning either `t.get()` or `t`
  depending on whether `t` is a  `std::reference_wrapper` or not.
- The header treedata.hh and the class `TreeData` are deprecated.
  Please use `TreeContainer` from dune-typetree instead.
- Function `derivative()` of `LocalAnalyticGridFunction` returns a bound
  local-function if the underlying value-local-function was bound before, i.e.,
  the bound-to element and geometry are passed to the derivative local-function.
- The concepts for `LocalFunction` and `LocalView` were extended with
  a `bound()` member function. It tells whether the object was already
  bound with a call to the `bind()` method or not.  Implementations of
  the two concepts outside dune-functions will have to implement this
  additional method.
- The derivative of a `DiscreteGlobalBasisFunction` is now implemented.

## Release 2.8

- `PreBasis` implementations are now required to provide a method
  `PreBasis::indices(node,iterator)` that replaces binding a `NodeIndexSet`
  to `node` and then calling `NodeIndexSet::indices(iterator)`.
  As a consequence `PreBasis::IndexSet` and `PreBasis::makeIndexSet`
  are no longer needed.

- The `RaviartThomasBasis` and `BrezziDouglasMariniBasis` now return
  Piola-transformed shape functions.  This is implemented by changing
  the return value of `tree().finiteElement()`: It is not an object
  of type `RaviartThomasLocalFiniteElement`
  or `BrezziDouglasMariniLocalFiniteElement` anymore. Rather, it is
  an object of a new type `GlobalValuedLocalFiniteElement`, which wraps
  other `LocalFiniteElement` implementations and applies a range-space
  transformation.  Domain-space transformations still have to be done
  by the calling code. The `GlobalValuedLocalFiniteElement` still
  implements the `LocalFiniteElement` interface of `dune-localfunctions`.

- The `RaviartThomasBasis` class now supports tetrahedral grids for `order=0`,
  quadrilateral grids for `order=2`, and hexahedral grids for `order=1`.

- The `RannacherTurekBasis` class now supports Crouzeix-Raviart elements.
  Grids containing simplices, cubes or both in 2d and 3d are supported now.

- The `dune-functions` module now contains an implementation of a
  Nedelec basis (for problems posed in H(curl)).  While the interface
  caters to different basis orders, grid dimensions and element types,
  only the first-order basis called "of the first kind" is implemented,
  and only for grids containing simplices, cubes or both in 2d and 3d.

- The `dune-functions` module now contains an implementation of a Hierarchical Lagrange
  basis for second order on simplex grids

- There is now an experimental implementation of a periodic basis in form
  of a class `PeriodicBasis`.  It is a meta basis, i.e., a basis that is
  parametrized with another basis (the host basis).  In a `PeriodicBasis`,
  global degrees of freedom of the host basis can be grouped into
  equivalence classes, which are then treated as single global degrees
  of freedom.  This allows, in particular, to implement periodic
  boundary conditions for discretizations without intersection integrals.

  The `PeriodicBasis` class can only be constructed by using the `periodic`
  method from the namespace `Dune::Functions::BasisFactory::Experimental`.
  It can change at any moment without much advance notice.  Use it at your
  own risk, and give us feedback!

- Imported the Python bindings from the 2.7 branch of dune-python and fixed remaining issues.
  Added a CI test that builds various global bases in 2d and 3d  and verifies the correct number of dofs.

- `interpolate` is now capable of interpolating vector-valued finite element functions correctly.
  The method of using scalar basis functions combined with vector-valued coefficients to mock a power basis is still supported.

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
