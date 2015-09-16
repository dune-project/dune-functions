add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        analyticgridviewfunction.hh
        discretescalarglobalbasisfunction.hh
        discreteglobalbasisfunction.hh
        gridfunction.hh
        gridfunction_imp.hh
        gridviewentityset.hh
        gridviewfunction.hh
        localderivativetraits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/gridfunctions)
