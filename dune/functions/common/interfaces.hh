// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_INTERFACES_HH
#define DUNE_FUNCTIONS_COMMON_INTERFACES_HH

#include <type_traits>


namespace Dune {
namespace Functions {



template<class Interface>
class PolymorphicType
{
public:

  /** \brief Destructor */
  virtual ~PolymorphicType()
  {}

  /** \brief Clones the object
   * clone needs to be redefined by an implementation class, with the
   * return type covariantly adapted. Remember to
   * delete the resulting pointer.
   */
  virtual Interface* clone() const = 0;

  virtual Interface* clone(void* buffer) const = 0;

  virtual Interface* move(void* buffer) = 0;
};



}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_EXTENDED_CLONEABLE_HH
