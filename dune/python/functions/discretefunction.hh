// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_PYTHON_FUNCTIONS_DISCRETEFUNCTION_HH
#define DUNE_PYTHON_FUNCTIONS_DISCRETEFUNCTION_HH

#include <dune/common/ftraits.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/python/common/pythonvector.hh>
#include <dune/python/grid/function.hh>
#include <dune/python/functions/hierarchicvectorwrapper.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {


    // HierarchicPythonVector
    // ----------------------

    template< class K >
    using HierarchicPythonVector = HierarchicVectorWrapper< PythonVector< K >, K >;



    namespace detail
    {

      // registerDiscreteFunctionConstructor
      // -----------------------------------

      template< class Basis, class K, class Range, class... options >
      inline static void registerDiscreteFunctionConstructor ( pybind11::class_< Dune::Functions::DiscreteGlobalBasisFunction< Basis, HierarchicPythonVector< K >, Dune::Functions::HierarchicNodeToRangeMap, Range >, options... > &cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;

        typedef HierarchicPythonVector< K > Vector;
        typedef Dune::Functions::HierarchicNodeToRangeMap NodeToRangeMap;
        typedef Dune::Functions::DiscreteGlobalBasisFunction< Basis, Vector, NodeToRangeMap, Range > DiscreteFunction;

        cls.def( pybind11::init( [] ( const Basis &basis, pybind11::buffer dofVector ) {
              auto nodeToRangeMapPtr =
                std::make_shared<const NodeToRangeMap>();
              std::shared_ptr<const Basis> basisPtr = Dune::wrap_or_move( basis );
              auto vectorPtr = std::make_shared< const Vector >( dofVector );
              return new DiscreteFunction( basisPtr,
                                           vectorPtr,
                                           nodeToRangeMapPtr );
            } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), "basis"_a, "dofVector"_a );
      }

      template< class DiscreteFunction, class... options >
      inline static void registerDiscreteFunctionConstructor ( pybind11::class_< DiscreteFunction, options... > &cls, PriorityTag< 0 > )
      {}

    } // namespace detail



    // registerDiscreteFunction
    // ------------------------

    template< class Basis, class Vector, class NodeToRangeMap, class Range, class... options >
    inline static void registerDiscreteFunction ( pybind11::module module, pybind11::class_< Dune::Functions::DiscreteGlobalBasisFunction< Basis, Vector, NodeToRangeMap, Range >, options... > &cls )
    {
      typedef Dune::Functions::DiscreteGlobalBasisFunction< Basis, Vector, NodeToRangeMap, Range > DiscreteFunction;

      registerGridFunction( module, cls );

      detail::registerDiscreteFunctionConstructor( cls, PriorityTag< 42 >() );

      cls.def_property_readonly( "basis", [] ( const DiscreteFunction &self ) -> const Basis & { return self.basis(); }, pybind11::keep_alive< 0, 1 >() );
      cls.def_property_readonly( "grid", [] ( const DiscreteFunction &self ) -> const typename Basis::GridView & { return self.basis().gridView(); }, pybind11::keep_alive< 0, 1 >() );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_FUNCTIONS_DISCRETEFUNCTION_HH
