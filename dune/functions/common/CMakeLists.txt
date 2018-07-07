add_subdirectory("test")

install(FILES
        callable.hh
        defaultderivativetraits.hh
        differentiablefunction.hh
        differentiablefunction_imp.hh
        differentiablefunctionfromcallables.hh
        functionfromcallable.hh
        functionconcepts.hh
        indexaccess.hh
        interfaces.hh
        lagrangelfecache.hh
        localfunction.hh
        localfunction_imp.hh
        optional.hh
        polymorphicsmallobject.hh
        reserveddeque.hh
        signature.hh
        staticforloop.hh
        staticgeometrytype.hh
        treedata.hh
        type_traits.hh
        typeerasure.hh
        utility.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
