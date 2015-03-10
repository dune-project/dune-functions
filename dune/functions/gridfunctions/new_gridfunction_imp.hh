// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_IMP_HH

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/interfaces.hh>



namespace Dune {
namespace Functions {
namespace Imp {



template<class Signature, class DerivativeInterface, class LocalFunctionInterface, class EntitySet>
class GridFunctionWrapperBase
{};

template<class Range, class Domain, class DerivativeInterface, class LocalFunctionInterface, class EntitySet>
class GridFunctionWrapperBase<Range(Domain), DerivativeInterface, LocalFunctionInterface, EntitySet> :
  public PolymorphicType<GridFunctionWrapperBase<Range(Domain), DerivativeInterface, LocalFunctionInterface, EntitySet> >
{
public:

  virtual Range operator() (const Domain& x) const = 0;

  virtual DerivativeInterface derivative() const = 0;

  virtual LocalFunctionInterface wrappedLocalFunction() const = 0;
};



template<class Signature, class DerivativeInterface, class LocalFunctionInterface, class EntitySet, class FImp>
class GridFunctionWrapper
{};

template<class Range, class Domain, class DerivativeInterface, class LocalFunctionInterface, class EntitySet, class FImp>
class GridFunctionWrapper< Range(Domain), DerivativeInterface, LocalFunctionInterface, EntitySet, FImp> :
  public GridFunctionWrapperBase<Range(Domain), DerivativeInterface, LocalFunctionInterface, EntitySet>
{
public:

  template<class F, disableCopyMove<GridFunctionWrapper, F> = 0>
  GridFunctionWrapper(F&& f) :
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

  virtual LocalFunctionInterface wrappedLocalFunction() const
  {
    return localFunction(f_);
  }

  virtual GridFunctionWrapper* clone() const
  {
    return new GridFunctionWrapper(*this);
  }

  virtual GridFunctionWrapper* clone(void* buffer) const
  {
    return new (buffer) GridFunctionWrapper(*this);
  }

  virtual GridFunctionWrapper* move(void* buffer)
  {
    return new (buffer) GridFunctionWrapper(std::move(*this));
  }

private:
  FImp f_;
};



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_IMP_HH
