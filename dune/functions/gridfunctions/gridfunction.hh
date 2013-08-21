// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDFUNCTION_HH

#include <memory>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/localfunction.hh>


namespace Dune {

namespace Functions {

/**
 * These classes describe the interface of EntitySets
 * that can be used as parameter for GridFunction.
 *
class EntitySet
{
  public:
    //! Type of Elements contained in this EntitySet
    typedef ... Element;

    //! Type of local coordinates with respect to the Element
    typedef Element::Geometry::LocalCoordinate LocalCoordinate;
    typedef Element::Geometry::GlobalCoordinate GlobalCoordinate;

    typedef Element value_type;

    //! Returns true if e is contained in the EntitySet
    bool contains(const Element& e) const;
};

class IterableEntitySet :
  public EntitySet
{
  public:

    //! A forward iterator
    typedef ... const_iterator;

    //! Same as const_iterator
    typedef const_iterator iterator;

    //! Number of Elements visited by an iterator
    size_t size() const;

    //! Create a begin iterator
    const_iterator begin() const;

    //! Create a end iterator
    const_iterator end() const;
};
*/


/** \brief Abstract base class for functions defined on a subset of grid entities.
 *
 * Being defined on a subset of entities means in particular that you can evaluate the function
 * in local coordinates of a given element contained in this subset. The given type EntitySet
 * is required to allow cheap copies and to implement the interface of EntitySet.
 *
 * \tparam ES Implementation of an EntitySet.
 * \tparam RT The type used for function values
 */
template<typename ES, typename RT>
class GridFunction
  : public DifferentiableFunction<typename ES::GlobalCoordinate,RT>
{

  typedef DifferentiableFunction<typename ES::GlobalCoordinate, RT> Base;

public:

  typedef ES EntitySet;
  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;
  typedef typename Base::DerivativeRange DerivativeRange;

  typedef GridFunction<ES,DerivativeRange> Derivative;

  typedef typename ES::LocalCoordinate LocalDomain;
  typedef typename ES::Element Element;

  typedef ::Dune::Functions::LocalFunction<GridFunction,Element> LocalFunction;

protected:

  typedef shared_ptr<LocalFunction> LocalFunctionBasePointer;

public:

  /** \brief Construction from a given EnitySet */
  GridFunction(const EntitySet& es)
    : es_(es)
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
  virtual Derivative* derivative() const = 0;

  /** \brief Const access to the grid view that the function is defined on */
  const EntitySet& entitySet() const
  {
    return es_;
  }

private:

  EntitySet es_;

};


} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDFUNCTION_HH
