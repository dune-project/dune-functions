// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH

#include <array>

#include <dune/common/hash.hh>

namespace Dune {
namespace Functions {



/**
 * \brief A multi index class with only one level
 *
 * This only adds a cast to size_type to std::array<size_type,1>.
 * Hence MultiIndices of type FlatMultiIndex can be used like
 * classic indices.
 */
template<class size_type>
class FlatMultiIndex :
  public std::array<size_type,1>
{
public:

  /**
   * \brief Forward constructor arguments to std::array
   */
  template<class... T>
  constexpr FlatMultiIndex(T&&... t) :
    std::array<size_type,1>(t...)
  {}

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

  operator const size_type& () const
  {
    return this->operator[](0);
  }

  operator size_type& ()
  {
    return this->operator[](0);
  }

  inline friend std::size_t hash_value(const FlatMultiIndex& arg)
  {
    return std::hash<size_type>()(arg);
  }

};



} // end namespace Functions
} // end namespace Dune

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename size_type),DUNE_HASH_TYPE(Dune::Functions::FlatMultiIndex<size_type>))

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH
