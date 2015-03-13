add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        bsplinebasis.hh
        gridviewfunctionspacebasis.hh
        interpolate.hh
        pq1nodalbasis.hh
        pq2nodalbasis.hh
        pqknodalbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
