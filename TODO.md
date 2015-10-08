# Open questions
<!-- vi: set ft=mkd ts=8 sw=2 et sts=2: -->

Feel free to put your wishes and questions into this file.
It would be nice to prefix your statements with your initials
or one more letter if needed.

## Development
* CaG: Should we discuss the further development using a bug tracker?
  And put everythig below into requests?


## Function space basis interface

* OS: The TaylorHoodIndexSet class currently contains two members

    enum { multiIndexMaxSize = 2 };
    size_type dimension() const;

  which were never officially discussed.  Christian added them because
  something like this was needed, but we may want to discuss these
  methods at the next meeting.

* Here is another point that needs future discussion:
  I just added the first implementation of a basis of a DG space.
  This one uses Lagrangian shape functions.  I expect other DG spaces using
  other shape functions to appear eventually.  Here is the issue:  Currently,
  the LagrangeDGBasis has the shape function type hard-wired into the code.
  Any implementation for a DG space using a different set of shape functions
  will differ very little from it.  Hence one may consider implementing
  a single DGBasis, and making the shape functions a template parameter.
  However, it is not clear what exactly this template parameter should be.
  It cannot be a LocalFiniteElement, because that would force us to pick
  a grid element type at compile time (you currently could not use the
  VirtualLocalFiniteElement here, because the code does not properly set it up).
  You cannot expect a factory here either, because then you cannot currently
  hard-wire the element type if that is what you want.
  I don't think this is a difficult issue, it just needs a bit of discussion.



## Function interface


### EntitySet interface
* CaG: It seems that the EntitySet concept and GridViewEntitySet are
  candidates for dune-grid.



### GridFunction concept
* CaG: Should GridFunction export the grid type?
* CaG: If yes, should GridFunction store/export a pointer to the grid?
* CaG: While LocalFunction is defined in a grid-agnostig way with
  local coordinates from a 'LocalContext'. In contrast to this
  the naming of GridFunction is linked to Grids:
  entitySet(), LocalContext=decltype(entitySet())::Element.
  There seems to be no need for this and we could
  simply rename this in a grid-agnostic way if someone
  comes up with better names than LocalizableFunction,
  localContextSet()



### GridViewFunction concept
* CaG: Should GridViewFunction default-implement the global evaluate() method?
  This would be as simple as

    typedef typename GridView::Grid::template Codim<0>::EntityPointer ElementPointer;
    Dune::HierarchicSearch<typename GridView::Grid, EntitySet> hsearch(this->gridView().grid(), this->entitySet());
    ElementPointer e = hsearch.findEntity(x);
    ::Dune::Functions::localFunction(*this, *e)->evaluate(e->geometry().local(x), y);

  In priciple this would also work for GridFunction if it knows the grid
  and if Element=Entity<0>.
  Done: No.