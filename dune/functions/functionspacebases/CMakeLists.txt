add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        bsplinebasis.hh
        concepts.hh
        defaultglobalbasis.hh
        defaultglobalindexset.hh
        defaultlocalindexset.hh
        defaultlocalview.hh
        flatmultiindex.hh
        flatvectorbackend.hh
        gridviewfunctionspacebasis.hh
        hierarchicvectorwrapper.hh
        interpolate.hh
        lagrangedgbasis.hh
        pq1nodalbasis.hh
        pq2nodalbasis.hh
        pqknodalbasis.hh
        nodes.hh
        sizeinfo.hh
        tupletreepath.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
