add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        callable.hh
        concept.hh
        differentiablefunction.hh
        functionfromcallable.hh
        localfunction.hh
        smallobject.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
