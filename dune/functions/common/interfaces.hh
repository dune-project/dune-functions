// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_INTERFACES_HH
#define DUNE_FUNCTIONS_COMMON_INTERFACES_HH

#include <type_traits>


namespace Dune {
namespace Functions {



/**
 * \brief Base class with polymorphic type boiler plate code
 *
 * \ingroup Utility
 *
 * By deriving your polymorphic type from this class you
 * enforce a default interface for common operations like
 * clone, clone to buffer, and move.
 */
template<class Interface>
class PolymorphicType
{
public:

  /** \brief Destructor */
  virtual ~PolymorphicType()
  {}

  /**
   * \brief Clones the object
   *
   * clone() needs to be redefined by an implementation class, with the
   * return type covariantly adapted. This will return a new copy of *this
   * via a pointer to newly allocated memory.
   * Remember to delete the resulting pointer.
   */
  virtual Interface* clone() const = 0;

  /**
   * \brief Clones the object into buffer
   *
   * clone(buffer) needs to be redefined by an implementation class,
   * with the return type covariantly adapted. This will return a copy
   * of *this created in the given buffer using placement-new with copy construction.
   * You must not delete the returned pointer since it points
   * to the given buffer (however with the proper type instead of void*).
   */
  virtual Interface* clone(void* buffer) const = 0;

  /**
   * \brief Move object into buffer
   *
   * move(buffer) needs to be redefined by an implementation class,
   * with the return type covariantly adapted. This will return a copy
   * of *this created in the given buffer using placement-new with move construction.
   * You must not delete the returned pointer since it points
   * to the given buffer (however with the proper type instead of void*).
   */
  virtual Interface* move(void* buffer) = 0;
};



}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_EXTENDED_CLONEABLE_HH
