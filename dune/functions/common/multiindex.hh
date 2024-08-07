// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_MULTIINDEX_HH
#define DUNE_FUNCTIONS_COMMON_MULTIINDEX_HH

#include <cstddef>
#include <array>
#include <iostream>

#include <dune/common/hash.hh>


namespace Dune::Functions {



/**
 * \brief A statically sized MultiIndex type
 *
 * This is basically a std::array<size_type,n>.
 */
template<class size_type, std::size_t n>
class StaticMultiIndex :
  public std::array<size_type, n>
{
public:
  static constexpr std::size_t size() { return n; }
  static constexpr std::size_t max_size() { return n; }

  inline friend std::size_t hash_value(const StaticMultiIndex& v) noexcept {
    return hash_range(v.begin(), v.end());
  }

};



/**
 * \brief A statically sized MultiIndex type
 *
 * This is basically a std::array<size_type,1>.
 *
 * This is the specialization for size==1 which
 * additionally provides const and mutable casts
 * to a reference to the single contained digit.
 */
template<class size_type>
class StaticMultiIndex<size_type,1> :
  public std::array<size_type, 1>
{
public:

  static constexpr std::size_t size() { return 1; }
  static constexpr std::size_t max_size() { return 1; }


  operator const size_type& () const {
    return (*this)[0];
  }

  inline friend std::size_t hash_value(const StaticMultiIndex& v) noexcept {
    return hash_range(v.begin(), v.end());
  }

  operator size_type& () {
    return (*this)[0];
  }

};



template<typename Stream, class size_type, std::size_t n>
inline Stream& operator<<(Stream& stream, const StaticMultiIndex<size_type,n>& c) {
  for (const auto& ci : c)
    stream << ci << "  ";
  return stream;
}



} // namespace Dune::Functions

template<class size_type, std::size_t n>
struct std::tuple_size< Dune::Functions::StaticMultiIndex<size_type,n> >
  : std::integral_constant<std::size_t, n> { };

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename size_type, std::size_t n),DUNE_HASH_TYPE(Dune::Functions::StaticMultiIndex<size_type,n>))

#endif // DUNE_FUNCTIONS_COMMON_MULTIINDEX_HH
