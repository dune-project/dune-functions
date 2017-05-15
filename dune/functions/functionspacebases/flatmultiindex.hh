// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH

#include <array>

#include <dune/common/hash.hh>

namespace Dune {
namespace Functions {



/**
 * \brief A multi-index class with only one level
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This class provides a multi-index interface in the sense
 * that it has operator[] access to individual interfaces.
 * However, since it only supports flat indices of exactly
 * one level, it also has a cast of the multi-index to
 * this index.
 * This is obtianed by deriving from std::array<size_type,1>
 * and adding this cast.
 * Hence multi-indices of type FlatMultiIndex can be used like
 * classic indices.
 */
template<class size_type>
class FlatMultiIndex :
  public std::array<size_type,1>
{
public:

  constexpr FlatMultiIndex() = default;

  /**
   * \brief Construct from initializer_list
   *
   * This is needed because std::array does not have
   * a constructor from initializer list. Instead
   * the list initialization of an std::array is
   * an aggregate initialization and hence not
   * visible in the derived class.
   */
  FlatMultiIndex(std::initializer_list<size_type> const &l) :
    std::array<size_type,1>{{*l.begin()}}
  {}

  /**
   * \brief Const cast of multi-index to first entry
   *
   * \returns The first (and only) index of this multi-index
   */
  operator const size_type& () const
  {
    return this->operator[](0);
  }

  /**
   * \brief Non-const cast of multi-index to first entry
   *
   * \returns The first (and only) index of this multi-index
   */
  operator size_type& ()
  {
    return this->operator[](0);
  }

  /**
   * \brief Compute hash value for FlatMultiIndex
   *
   * \ingroup FunctionSpaceBasesUtilities
   */
  inline friend std::size_t hash_value(const FlatMultiIndex& arg) noexcept
  {
    return std::hash<size_type>()(arg);
  }

};



} // end namespace Functions
} // end namespace Dune

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename size_type),DUNE_HASH_TYPE(Dune::Functions::FlatMultiIndex<size_type>))

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH
