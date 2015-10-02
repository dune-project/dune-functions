add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        callable.hh
        concept.hh
        defaultderivativetraits.hh
        differentiablefunction.hh
        differentiablefunction_imp.hh
        differentiablefunctionfromcallables.hh
        functionfromcallable.hh
        functionconcepts.hh
        interfaces.hh
        localfunction.hh
        localfunction_imp.hh
        optional.hh
        polymorphicsmallobject.hh
        reserveddeque.hh
        signature.hh
        staticforloop.hh
        treedata.hh
        type_traits.hh
        typeerasure.hh
        typelist.hh
        utility.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
