add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        callable.hh
        concept.hh
        defaultderivativetraits.hh
        differentiablefunction.hh
        differentiablefunction_imp.hh
        functionfromcallable.hh
        interfaces.hh
        localfunction.hh
        localfunction_imp.hh
        signature.hh
        smallobject.hh
        type_traits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
