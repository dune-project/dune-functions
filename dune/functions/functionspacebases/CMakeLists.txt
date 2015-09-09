add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        bsplinebasis.hh
        defaultlocalindexset.hh
        defaultglobalindexset.hh
        gridviewfunctionspacebasis.hh
        interpolate.hh
        lagrangedgbasis.hh
        pq1nodalbasis.hh
        pq2nodalbasis.hh
        pqknodalbasis.hh
        nodes.hh
        tupletreepath.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
