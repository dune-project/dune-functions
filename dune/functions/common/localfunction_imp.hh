// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_LOCALFUNCTION_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_COMMON_LOCALFUNCTION_FUNCTION_IMP_HH

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/interfaces.hh>



namespace Dune {
namespace Functions {
namespace Imp {



template<class Signature, class DerivativeInterface, class LocalContext>
class LocalFunctionWrapperBase
{};

template<class Range, class Domain, class DerivativeInterface, class LocalContext>
class LocalFunctionWrapperBase<Range(Domain), DerivativeInterface, LocalContext> :
  public PolymorphicType<LocalFunctionWrapperBase<Range(Domain), DerivativeInterface, LocalContext> >
{
public:

  virtual Range operator() (const Domain& x) const = 0;

  virtual DerivativeInterface derivative() const = 0;

  virtual void bind(const LocalContext&) = 0;

  virtual void unbind() = 0;

  virtual const LocalContext& localContext() const = 0;
};



template<class Signature, class DerivativeInterface, class LocalContext, class FImp>
class LocalFunctionWrapper
{};

template<class Range, class Domain, class DerivativeInterface, class LocalContext, class FImp>
class LocalFunctionWrapper< Range(Domain), DerivativeInterface, LocalContext, FImp> :
  public LocalFunctionWrapperBase<Range(Domain), DerivativeInterface, LocalContext>
{
public:

  template<class F, disableCopyMove<LocalFunctionWrapper, F> = 0>
  LocalFunctionWrapper(F&& f) :
    f_(std::forward<F>(f))
  {}

  virtual Range operator() (const Domain& x) const
  {
    return f_(x);
  }

  virtual DerivativeInterface derivative() const
  {
    return derivativeIfImplemented<DerivativeInterface, FImp>(f_);
  }

  virtual void bind(const LocalContext& context)
  {
    return f_.bind(context);
  }

  virtual void unbind()
  {
    return f_.unbind();
  }

  virtual const LocalContext& localContext() const
  {
    return f_.localContext();
  }

  virtual LocalFunctionWrapper* clone() const
  {
    return new LocalFunctionWrapper(*this);
  }

  virtual LocalFunctionWrapper* clone(void* buffer) const
  {
    return new (buffer) LocalFunctionWrapper(*this);
  }

  virtual LocalFunctionWrapper* move(void* buffer)
  {
    return new (buffer) LocalFunctionWrapper(std::move(*this));
  }

private:
  FImp f_;
};



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
