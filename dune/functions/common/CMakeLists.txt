add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        callable.hh
        differentiablefunction.hh
        functionfromcallable.hh
        localfunction.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
