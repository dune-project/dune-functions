// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFPREBASISMIXIN_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFPREBASISMIXIN_HH

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace Dune::Functions {

/**
 * \brief A generic MixIn class for PreBasis
 *
 * Extends the interface of a `Derived` class by common
 * constants and the `size()` methods. A requirement is
 * that `Derived` implements the method `dimension()`
 * returning the total number of basis functions the
 * pre-basis represents.
 *
 * This mixin class can be used for all pre bases that
 * are on the leaf of the basis-tree. These pre-bases
 * are supposed to have a size only for empty size-prefixes
 * and this size is the same as the given dimension.
 *
 * \tparam Derived  The actual implementation of a pre-basis
 */
template<class Derived>
class LeafPreBasisMixin
{
public:
  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Maximal length of global multi-indices
  static constexpr size_type maxMultiIndexSize = 1;

  //! Minimal length of global multi-indices
  static constexpr size_type minMultiIndexSize = 1;

  //! Size required temporarily when constructing global multi-indices
  static constexpr size_type multiIndexBufferSize = 1;

  //! Return number of possible values for next position in multi index
  template<class SizePrefix,
    decltype(std::declval<SizePrefix>().size(), bool{}) = true>
  size_type size(const SizePrefix& prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? derived().dimension() : 0;
  }

  //! Get the total dimension of the space spanned by this basis
  size_type size() const
  {
    return derived().dimension();
  }

private:
  const Derived& derived() const
  {
    return static_cast<const Derived&>(*this);
  }
};


} // end namespace Dune::Functions


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFPREBASISMIXIN_HH
