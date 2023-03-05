add_subdirectory("test")

install(FILES
        defaultderivativetraits.hh
        differentiablefunction.hh
        differentiablefunction_imp.hh
        differentiablefunctionfromcallables.hh
        functionconcepts.hh
        indexaccess.hh
        interfaces.hh
        localfunction.hh
        localfunction_imp.hh
        multiindex.hh
        overflowarray.hh
        polymorphicsmallobject.hh
        reserveddeque.hh
        signature.hh
        staticforloop.hh
        type_traits.hh
        typeerasure.hh
        utility.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
