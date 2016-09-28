// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_POLYMORPHICSMALLOBJECT_HH
#define DUNE_FUNCTIONS_COMMON_POLYMORPHICSMALLOBJECT_HH

#include <utility>

namespace Dune {
namespace Functions {


/**
 * \brief A wrapper providing small object optimization with polymorphic types
 *
 * \ingroup Utility
 * \ingroup TypeErasure
 *
 * \tparam Base Base class type of wrapped objects
 * \tparam bufferSize Size of small object buffer
 *
 * This class encapsulates small object optimization for polymorphic types.
 * The type of objects passed to the constructor must be derived from the
 * base class type Base.
 *
 * If the size of the derived type fits into the static buffer, then the
 * wrapped object is stored there, otherwise it is allocated dynamically.
 *
 * In order to make the copy constructor work for polymorphic types,
 * Base must provide clone() and clone(void*). The former should return
 * a pointer to a dynamically allocated clone, while the latter
 * should call the appropriate placement-new with the passed pointer.
 *
 * Similarly the polymorphic type has to implement move(void*).
 * This should call placement-new and can std::move all the
 * data but leave the object in a valid and probably unusable state.
 */
template<class Base, size_t bufferSize>
class PolymorphicSmallObject
{
public:

  //! Default constructor
  PolymorphicSmallObject() :
    p_(nullptr)
  {}

  /**
   * \brief Construct from object
   *
   * \tparam Derived Type of object to be stored, should be derived from Base
   * \param derived Object to be stored
   */
  template<class Derived>
  PolymorphicSmallObject(Derived&& derived)
  {
    if (sizeof(Derived)<bufferSize)
      p_ = new (buffer_) Derived(std::forward<Derived>(derived));
    else
      p_ = new Derived(std::forward<Derived>(derived));
  }

  //! Move constructor from other PolymorphicSmallObject
  PolymorphicSmallObject(PolymorphicSmallObject&& other)
  {
    moveToWrappedObject(std::move(other));
  }

  //! Copy constructor from other PolymorphicSmallObject
  PolymorphicSmallObject(const PolymorphicSmallObject& other)
  {
    copyToWrappedObject(other);
  }

  //! Destructor
  ~PolymorphicSmallObject()
  {
    destroyWrappedObject();
  }

  //! Copy assignment from other PolymorphicSmallObject
  PolymorphicSmallObject& operator=(const PolymorphicSmallObject& other)
  {
    destroyWrappedObject();
    copyToWrappedObject(other);
    return *this;
  }

  //! Move assignment from other PolymorphicSmallObject
  PolymorphicSmallObject& operator=(PolymorphicSmallObject&& other)
  {
    destroyWrappedObject();
    moveToWrappedObject(std::move(other));
    return *this;
  }

  //! Check if *this is not empty
  explicit operator bool() const
  {
    return p_;
  }

  //! Check if object is stored in internal stack buffer
  bool bufferUsed() const
  {
    return ((void*) (p_) == (void*)(&buffer_));
  }

  //! Obtain reference to stored object
  const Base& get() const
  {
    return *p_;
  }

  //! Obtain mutable reference to stored object
  Base& get()
  {
    return *p_;
  }

private:

  void destroyWrappedObject()
  {
    if (operator bool())
    {
      if (bufferUsed())
        p_->~Base();
      else
        delete p_;
    }
  }

  void moveToWrappedObject(PolymorphicSmallObject&& other)
  {
    if (other.bufferUsed())
      p_ = other.p_->move(buffer_);
    else
    {
      // We don't need to check for &other_!=this, because you can't
      // have an rvalue to *this and call it's assignment/constructor
      // at the same time. (Despite trying to shot yourself in the foot
      // with std::move explicitly.)

      // Take ownership of allocated object
      p_ = other.p_;

      // Leave pointer in a clear state to avoid double freeing it.
      other.p_ = 0;
    }
  }

  void copyToWrappedObject(const PolymorphicSmallObject& other)
  {
    if (&other!=this)
    {
      if (other.bufferUsed())
        p_ = other.p_->clone(buffer_);
      else
        p_ = other.p_->clone();
    }
  }

  alignas(Base) char buffer_[bufferSize];
  Base* p_;
};


} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_POLYMORPHICSMALLOBJECT_HH
