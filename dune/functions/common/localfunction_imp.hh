// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_LOCALFUNCTION_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_COMMON_LOCALFUNCTION_FUNCTION_IMP_HH

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>



namespace Dune {
namespace Functions {
namespace Imp {

// Interface of type erasure wrapper
//
// Notice that the basic interface of polymorphic classes (destructor, clone, ...)
// will be added by the type erasure foundation classes.
template<class Signature, class DerivativeInterface, class LocalContext>
class LocalFunctionWrapperInterface :
  public DifferentiableFunctionWrapperInterface<Signature, DerivativeInterface>
{
public:
  virtual void bind(const LocalContext&) = 0;

  virtual void unbind() = 0;

  virtual bool bound() const = 0;

  virtual const LocalContext& localContext() const = 0;
};


// Implementation of type erasure wrapper
template<class Signature, class DerivativeInterface, class LocalContext, class B>
class LocalFunctionWrapperImplementation :
  public DifferentiableFunctionWrapperImplementation<Signature, DerivativeInterface, B>
{
  using Base = DifferentiableFunctionWrapperImplementation<Signature, DerivativeInterface, B>;
public:
  using Base::Base;

  virtual void bind(const LocalContext& context)
  {
    this->get().bind(context);
  }

  virtual void unbind()
  {
    this->get().unbind();
  }

  virtual bool bound() const
  {
    return this->get().bound();
  }

  virtual const LocalContext& localContext() const
  {
    return this->get().localContext();
  }
};

}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
