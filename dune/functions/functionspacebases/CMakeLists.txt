add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        bsplinebasis.hh
        gridviewfunctionspacebasis.hh
        pq1nodalbasis.hh
        pq2nodalbasis.hh
        pq3nodalbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
