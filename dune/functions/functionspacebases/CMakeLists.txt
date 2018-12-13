add_subdirectory("test")

install(FILES
        basistags.hh
        boundarydofs.hh
        brezzidouglasmarinibasis.hh
        bsplinebasis.hh
        compositebasis.hh
        concepts.hh
        defaultglobalbasis.hh
        defaultlocalview.hh
        defaultnodetorangemap.hh
        flatmultiindex.hh
        flatvectorview.hh
        hierarchicnodetorangemap.hh
        hierarchicvectorwrapper.hh
        interpolate.hh
        lagrangebasis.hh
        lagrangedgbasis.hh
        powerbasis.hh
        rannacherturekbasis.hh
        nodes.hh
        sizeinfo.hh
        subentitydofs.hh
        subspacebasis.hh
        subspacelocalview.hh
        taylorhoodbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
