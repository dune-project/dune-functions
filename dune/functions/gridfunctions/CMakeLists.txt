add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        discretescalarglobalbasisfunction.hh
        gridfunction.hh
        gridfunction_imp.hh
        gridviewentityset.hh
        gridviewfunction.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/gridfunctions)
