<!-- vi: set ft=mkd ts=8 sw=2 et sts=2: -->
# The Dune-Function module

## Scope of the module

The dune-functions module provides an abstraction layer for global finite
element functions. Its two main concepts are functions implemented as callable
objects, and bases of finite element spaces.

### Functions

dune-functions provides an interface to "functions" in the mathematical sense,
in particular to finite element functions defined on a grid, but going far
beyond that.

The interface revolves around the concept of a "callable". It encompasses any
type of C++ object that can be evaluated with `operator()`, like free functions,
function objects, and even C++11 lambdas. Dynamic polymorphism is realized
using type erasure and the `std::function` class, which does not sacrifice
efficiency in purely static code.

dune-functions extends the "callable" concept into several directions. First,
it allows for differentiable functions. Such functions can hand out their
derivative as new function objects. Second, for functions defined piecewisely
on a finite element grid, the concept of local function is introduced. Local
functions can be bound to grid elements. All further evaluations of a function
bound to an element are in local coordinates of that element. This approach
allows to avoid overhead when there are many evaluations on a single element.

For more details refer to the @ref Functions section.


### Function space bases

The second part of dune-functions provides a well-defined interface to bases of
finite element function spaces. For this interface, a finite element basis is a
set of functions with a prescribed ordering, and a way to index them. The core
functionality has three parts:

1.  For a given grid element, obtain the restrictions of all basis functions to
    this element, except for those functions where the restriction is zero. In
    other words: get the shape functions for the element.
2.  Get a local numbering for these shape functions. This is needed to index the element stiffness matrix.
3.  Get a global numbering for the shape functions. This is needed to index the global stiffness matrix.

While local numbers are always integers, global numbers can be multi-indices,
if appropriate.

A central feature is that finite element bases for vector-valued and mixed
spaced can be constructed by tensor multiplication of simpler bases. The
resulting expressions can be interpreted as tree structures. For example, the
tree for the three-dimensional Taylor-Hood basis is shown above. This tree
structure is directly exposed in the dune-functions interface. An easy
mechanism allows to construct new spaces.

For more details refer to the @ref FunctionSpaceBases section.

## Documentation

### Class documentation
The module contains a class documentation which can be build using [doxygen].
After the module has been build, you can build the documentation using
`make doc`
Additionally the pre-build doxygen documentation for the _master_ and
release branches is also hosted on the [documentation section][dune docs]
of the dune-website.

### Manual
There are two documents describing the concepts and usage of dune functions.
The interface of function is described in the article

    C. Engwer, C. Gräser, S. Müthing, and O. Sander.
    The interface for functions in the dune-functions module.
    Archive of Numerical Software, 5(1):95--109, 2017.

This is freely available
via the [website of the journal][functions paper] and
as [arXiv:1512.06136][functions paper arxiv] preprint.
The interface of the function space bases is described in the article

    C. Engwer, C. Gräser, S. Müthing, and O. Sander.
    Function space bases in the dune-functions module.
    Preprint, arxiv:1806.09545, 2018.

This is freely available
as [arXiv:1806.09545][bases paper arxiv] preprint.
Both are also contained in the module. Like the class documentation
this is build on `make doc`.

### Examples
Several example applications demonstrate how to use the module. These
example applications are contained in the `examples/` directory and
build when building the module. The `stokes-taylorhood` example is
described in detail in the manual (see above).


## Using dune-functions and licensing
The module is licensed by different variants of the GPL licence.
Please have a look at the `COPYING` file for more information
and a list of all contributors. When using dune-functions
**please make sure to cite the publications on the
[functions interface][functions paper] and the
[bases interface][bases paper]** listed above.



## Building dune-functions

### Dependencies
Dune-functions depends on the dune [core modules][core]
and the [dune-typetree module][typetree]. All of them are available using git:

* https://gitlab.dune-project.org/core/dune-common
* https://gitlab.dune-project.org/core/dune-geometry
* https://gitlab.dune-project.org/core/dune-grid
* https://gitlab.dune-project.org/core/dune-istl
* https://gitlab.dune-project.org/core/dune-localfunctions
* https://gitlab.dune-project.org/staging/dune-typetree

The versioning of dune-functions follows the scheme used in the core modules.
I.e. version x.y of dune-functions will depend on version x.y of the core modules
and dune-typetree. Analogously, the _master_ branch will depend on the
_master_ branch of these modules.

Unless explicitly stated otherwise for a specific version,
dune-functions supports/requires the same build tools (compilers, cmake)
as the corresponding version of the core modules.

### Building the module
Dune-functions integrates into the cmake-based dune build system.
Hence it can be build (like any other module) using the `dunecontrol` script
provided by the core modules. For details on how to use this build system
and how to specify build options have a look at the documentation in the
dune-common module.


[core]: https://dune-project.org/groups/core
[typetree]: https://gitlab.dune-project.org/staging/dune-typetree
[dune docs]: https://dune-project.org/doxygen
[functions paper arxiv]: https://arxiv.org/abs/1512.06136
[functions paper]: http://journals.ub.uni-heidelberg.de/index.php/ans/article/view/27683
[bases paper arxiv]: https://arxiv.org/abs/1806.09545
[bases paper]: https://arxiv.org/abs/1806.09545
[doxygen]: http://www.stack.nl/~dimitri/doxygen/
