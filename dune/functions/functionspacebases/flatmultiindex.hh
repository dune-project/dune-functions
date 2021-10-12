// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH

#include <array>

#include <dune/functions/common/multiindex.hh>

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
using FlatMultiIndex = StaticMultiIndex<size_type, 1>;



} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FLATMULTIINDEX_HH
