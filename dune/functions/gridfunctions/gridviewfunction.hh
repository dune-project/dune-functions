// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH

#include <memory>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/localfunction.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>


namespace Dune {

namespace Functions {

/** \brief Abstract base class for functions defined on a GridView
 *
 * Being defined on a grid view means in particular that you can evaluate the function
 * in local coordinates of a given element of the grid view.
 *
 * \tparam GV The GridView that the function is defined on
 * \tparam RT The type used for function values
 */
template<typename GV, typename RT>
class GridViewFunction
  : public GridFunction<GridViewEntitySet<GV, 0>, RT>
{

  typedef GridFunction<GridViewEntitySet<GV, 0>, RT> Base;

public:

  typedef GV GridView;
  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;
  typedef typename Base::DerivativeRange DerivativeRange;
  typedef typename Base::LocalDomain LocalDomain;
  typedef typename Base::Element Element;
  typedef typename Base::EntitySet EntitySet;

  typedef GridViewFunction<GV,DerivativeRange> Derivative;
  typedef typename Base::DerivativeBasePointer DerivativeBasePointer;

  typedef typename ::Dune::Functions::LocalFunction<Base,Element> LocalFunction;

protected:

  typedef shared_ptr<LocalFunction> LocalFunctionBasePointer;

public:

  /** \brief Construction from a given grid view */
  GridViewFunction(const GridView& gv)
    : Base(EntitySet(gv))
  {}

  /** \brief Access to the function on a single element, in coordinates of that element
   *
   * To evaluate the function on a single element you have to get a local function for
   * this element.  You can do this by calling this function and then binding the
   * local function to a given element.  Then the local function can be evaluated
   * at given points.
   *
   * Rationale: if you want to evaluate the function at many points in the same element
   * this approach is more efficient.
   */
  virtual LocalFunctionBasePointer localFunction() const = 0;

  /** \brief Access to the derivative function
   *
   * We pretend that the function is differentiable everywhere, even though this will
   * usually only be true in the interiors of the elements.
   */
  virtual typename Base::DerivativeBasePointer derivative() const = 0;

  /** \brief Const access to the grid view that the function is defined on */
  const GridView& gridView() const
  {
    return this->entitySet().gridView();
  }

private:

};


} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWFUNCTION_HH
