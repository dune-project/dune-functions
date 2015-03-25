// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_IMP_HH

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/interfaces.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>



namespace Dune {
namespace Functions {
namespace Imp {

/**
 * A concept describing types that have a localFunction() method found by ADL
 */
struct HasFreeLocalFunction
{
  template<class F>
  auto require(F&& f) -> decltype(
    localFunction(f)
  );
};



// Interface of type erasure wrapper
//
// Notice that the basic interface of polymorphic classes (destructor, clone, ...)
// will be added by the type erasure foundation classes.
template<class Signature, class DerivativeInterface, class LocalFunctionInterface, class EntitySet>
class GridFunctionWrapperInterface :
  public DifferentiableFunctionWrapperInterface<Signature, DerivativeInterface>
{
public:
  virtual LocalFunctionInterface wrappedLocalFunction() const = 0;

  virtual const EntitySet& wrappedEntitySet() const = 0;
};


// Implementation of type erasure wrapper
template<class Signature, class DerivativeInterface, class LocalFunctionInterface, class EntitySet, class B>
class GridFunctionWrapperImplementation :
  public DifferentiableFunctionWrapperImplementation<Signature, DerivativeInterface, B>
{
  using Base = DifferentiableFunctionWrapperImplementation<Signature, DerivativeInterface, B>;
public:
  using Base::Base;

  virtual LocalFunctionInterface wrappedLocalFunction() const
  {
    return localFunction(this->get());
  }

  virtual const EntitySet& wrappedEntitySet() const
  {
    return this->get().entitySet();
  }
};



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_IMP_HH
