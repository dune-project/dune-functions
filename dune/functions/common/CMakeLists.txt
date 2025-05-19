# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

add_subdirectory("test")

install(FILES
        defaultderivativetraits.hh
        densevectorview.hh
        differentiablefunction.hh
        differentiablefunction_imp.hh
        differentiablefunctionfromcallables.hh
        functionconcepts.hh
        geometryinancestor.hh
        indexaccess.hh
        interfaces.hh
        localfunction.hh
        localfunction_imp.hh
        mapperutilities.hh
        multiindex.hh
        overflowarray.hh
        polymorphicsmallobject.hh
        reserveddeque.hh
        signature.hh
        staticforloop.hh
        subdomain.hh
        type_traits.hh
        typeerasure.hh
        utility.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
