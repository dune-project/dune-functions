// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TYPEERASURE_HH
#define DUNE_FUNCTIONS_COMMON_TYPEERASURE_HH

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/interfaces.hh>

namespace Dune {
namespace Functions {
namespace Imp {



template<class Interface>
class TypeErasureBase :
  public Interface,
  public PolymorphicType<TypeErasureBase<Interface>>
{};



template<class Interface, class T>
class TypeErasureImp :
  public TypeErasureBase<Interface>
{
public:
  template<class TT, disableCopyMove<TypeErasureImp, TT> = 0>
  TypeErasureImp(TT&& t) :
    wrapped_(std::forward<TT>(t))
  {}

protected:
  using Wrapped = T;
  Wrapped wrapped_;
};



template<class Interface, template<class> class Implementation, class T>
class TypeErasureDerived :
  public Implementation<TypeErasureImp<Interface, T> >
{
public:

  template<class TT, disableCopyMove<TypeErasureDerived, T> = 0>
  TypeErasureDerived(TT&& t) :
    Implementation<TypeErasureImp<Interface, T> >(std::forward<TT>(t))
  {}

  virtual TypeErasureDerived* clone() const
  {
    return new TypeErasureDerived(*this);
  }

  virtual TypeErasureDerived* clone(void* buffer) const
  {
    return new (buffer) TypeErasureDerived(*this);
  }

  virtual TypeErasureDerived* move(void* buffer)
  {
    return new (buffer) TypeErasureDerived(std::move(*this));
  }
};



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_COMMON_TYPEERASURE_HH
