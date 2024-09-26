# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

add_subdirectory("test")

install(FILES
        basistags.hh
        boundarydofs.hh
        brezzidouglasmarinibasis.hh
        bsplinebasis.hh
        compositebasis.hh
        concepts.hh
        containerdescriptors.hh
        defaultglobalbasis.hh
        defaultlocalview.hh
        defaultnodetorangemap.hh
        dynamicpowerbasis.hh
        flatmultiindex.hh
        flatvectorview.hh
        globalvaluedlocalfiniteelement.hh
        hermitebasis.hh
        hierarchicallagrangebasis.hh
        hierarchicnodetorangemap.hh
        hierarchicvectorwrapper.hh
        interpolate.hh
        lagrangebasis.hh
        lagrangedgbasis.hh
        leafprebasismappermixin.hh
        leafprebasismixin.hh
        lfeprebasismixin.hh
        nedelecbasis.hh
        periodicbasis.hh
        powerbasis.hh
        rannacherturekbasis.hh
        raviartthomasbasis.hh
        refinedlagrangebasis.hh
        nodes.hh
        sizeinfo.hh
        subentitydofs.hh
        subspacebasis.hh
        subspacelocalview.hh
        taylorhoodbasis.hh
        transformedfiniteelementmixin.hh
        transformedindexbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
