// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZEINFO_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZEINFO_HH

#warning The header dune/functions/functionspacebases/sizeinfo.hh is deprecated and will be removed after release 2.10.

#include <array>

namespace Dune {
namespace Functions {



/**
 * \brief A class encapsulating size information
 *
 * This class encapsulates size information of a basis.
 *
 * We may want to return this or a similar class when
 * calling basis.size() without arguments. Until we
 * decided on this, we can use sizeInfo(basis) and
 * use the result as size provider for a VectorBackend.
 *
 * \deprecated This class is deprecated.
 */
template<class B>
class
[[deprecated("SizeInfo is deprecated and will be removed after release 2.10.")]]
SizeInfo
{
public:
    using Basis = B;
    using size_type = typename Basis::size_type;
    using SizePrefix = typename Basis::SizePrefix;

    /**
     * \brief Construct from basis
     */
    SizeInfo(const Basis& basis) :
        basis_(&basis)
    {}

    /**
     * \brief Return number possible values for next position in multi index
     */
    size_type operator()(const SizePrefix& prefix) const
    {
        return basis_->size(prefix);
    }

    /**
     * \brief Return number possible values for next position in multi index
     *
     * This shall vanish. It's just here such that this can be used
     * as size provider n place of the basis.
     */
    size_type size(const SizePrefix& prefix) const
    {
        return basis_->size(prefix);
    }

    operator size_type () const
    {
        return basis_->dimension();
    }

protected:

    const Basis* basis_;
};


/**
 * \deprecated This function is deprecated.
 */
template<class Basis>
SizeInfo<Basis> sizeInfo(const Basis& basis)
{
    return SizeInfo<Basis>(basis);
}



} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZEINFO_HH
