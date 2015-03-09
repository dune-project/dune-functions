#add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        polynomial.hh
        trigonometricfunction.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/analyticfunctions)
