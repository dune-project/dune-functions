add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        gridviewfunctionspacebasis.hh
        pq1nodalbasis.hh
        pq2nodalbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
