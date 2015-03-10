add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        gridfunction.hh
        gridviewentityset.hh
        gridviewfunction.hh
        new_gridfunction.hh
        new_gridfunction_imp.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/gridfunctions)
