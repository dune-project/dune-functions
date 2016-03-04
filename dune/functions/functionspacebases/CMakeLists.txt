add_subdirectory("test")

install(FILES
        basistags.hh
        bsplinebasis.hh
        compositebasis.hh
        concepts.hh
        defaultglobalbasis.hh
        defaultglobalindexset.hh
        defaultlocalindexset.hh
        defaultlocalview.hh
        defaultnodetorangemap.hh
        flatmultiindex.hh
        flatvectorbackend.hh
        gridviewfunctionspacebasis.hh
        hierarchicvectorwrapper.hh
        interpolate.hh
        lagrangedgbasis.hh
        powerbasis.hh
        pq1nodalbasis.hh
        pqknodalbasis.hh
        nodes.hh
        sizeinfo.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
