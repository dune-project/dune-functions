<!-- vi: set ft=mkd ts=8 sw=2 et sts=2: -->
# The Dune-Function module


The dune-functions module provides an abstraction layer for global finite
element functions. Its two main concepts are functions implemented as callable
objects, and bases of finite element spaces.

## Functions

dune-functions provides an interface to "functions" in the mathematical sense,
in particular to finite element functions defined on a grid, but going far
beyond that.

The interface revolves around the concept of a "callable". It encompasses any
type of C++ object that can be evaluated with operator(), like free functions,
function objects, and even C++11 lambdas. Dynamic polymorphism is realized
using type erasure and the std::function class, which does not sacrifice
efficiency in purely static code.

dune-functions extends the "callable" concept into several directions. First,
it allows for differentiable functions. Such functions can hand out their
derivative as new function objects. Second, for functions defined piecewisely
on a finite element grid, the concept of local function is introduced. Local
functions can be bound to grid elements. All further evaluations of a function
bound to an element are in local coordinates of that element. This approach
allows to avoid overhead when there are many evaluations on a single element.

For more details refer to the @ref Functions section.


## Function space bases

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

