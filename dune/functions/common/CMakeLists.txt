add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        callable.hh
        concept.hh
        defaultderivativetraits.hh
        differentiablefunction.hh
        functionfromcallable.hh
        localfunction.hh
        interfaces.hh
        new_differentiablefunction.hh
        new_differentiablefunction_imp.hh
        new_localfunction.hh
        new_localfunction_imp.hh
        smallobject.hh
        type_traits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
