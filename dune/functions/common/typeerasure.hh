// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TYPEERASURE_HH
#define DUNE_FUNCTIONS_COMMON_TYPEERASURE_HH

#include <typeinfo>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/interfaces.hh>
#include <dune/functions/common/polymorphicsmallobject.hh>

namespace Dune {
namespace Functions {
namespace Imp {



/**
 * \brief The internal wrapper interface for type erasure
 *
 * This class adds some foundation interfaces needed
 * for proper dynamic polymorphism and type erasure.
 *
 * The actual application interface has to be provided
 * by the user.
 *
 * \tparam Interface Class defininig the internal abstract virtual interface
 */
template<class Interface>
class TypeErasureBase :
  public Interface,
  public PolymorphicType<TypeErasureBase<Interface>>
{
public:
  virtual const std::type_info& target_type() const = 0;
};



/**
 * \brief Base implementation of the internal wrapper interface
 *
 * This class derives from the foundation interfaces
 * and the user defined interfaces provided by the interface
 * parameter. It will store any suitable type T to do
 * the type erasure.
 *
 * The implementation of the foundation and user interface
 * is provided by classed derived of this one.
 *
 * \tparam Interface Class defininig the internal abstract virtual interface
 * \tparam T A type modelleding the desired interface
 */
template<class Interface, class T>
class TypeErasureImp :
  public TypeErasureBase<Interface>
{
public:
  template<class TT, disableCopyMove<TypeErasureImp, TT> = 0>
  TypeErasureImp(TT&& t) :
    wrapped_(std::forward<TT>(t))
  {}

  T& get()
  {
    return wrapped_;
  }

  const T& get() const
  {
    return wrapped_;
  }

protected:
  using Wrapped = T;
  Wrapped wrapped_;
};



/**
 * \brief Implementation of the internal wrapper interface
 *
 * This class implements the foundation and user interfaces
 * of the internal type erasure wrapper.
 *
 * The foundation interface of TypeErasureBase is directly
 * implemented here whereas the user interface is implemented
 * by deriving from the user-provides Implementation template.
 *
 * The Implementation is a template taking one class template
 * parameter. It should directly or indirectly derive from this
 * class and inherit its constructors.
 * In order to forward the implemented methods to the erased
 * type it can use the wrapper_ member of this base class being
 * of this type.
 *
 * \tparam Interface Class defininig the internal abstract virtual interface
 * \tparam Implementation Class defininig implemention the abstract methods of Interface
 * \tparam T A type modelleding the desired interface
 */
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

  virtual const std::type_info& target_type() const
  {
    return typeid(T);
  }
};

} // namespace Dune::Functions::Imp

template<class Interface, template<class> class Implementation, size_t bufferSize = 56>
class TypeErasure
{
public:

  template<class T, disableCopyMove<TypeErasure, T> = 0 >
  TypeErasure(T&& t) :
    wrapped_(Imp::TypeErasureDerived<Interface, Implementation, typename std::decay<T>::type>(std::forward<T>(t)))
  {}

  TypeErasure() = default;

  Interface& asInterface()
  {
    return wrapped_.get();
  }

  const Interface& asInterface() const
  {
    return wrapped_.get();
  }

  virtual const std::type_info& target_type() const
  {
    return wrapped_.get().target_type();
  }

protected:
  PolymorphicSmallObject<Imp::TypeErasureBase<Interface>, bufferSize > wrapped_;
};


}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_TYPEERASURE_HH
