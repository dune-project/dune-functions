#ifndef DUNE_PYTHON_FUNCTIONS_SUBSPACEBASIS_HH
#define DUNE_PYTHON_FUNCTIONS_SUBSPACEBASIS_HH

#include <type_traits>

#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dune/python/functions/globalbasis.hh>

namespace Dune
{

  namespace Python
  {

    template< class SubspaceBasis, class... options >
    DUNE_EXPORT void registerSubspaceBasis ( pybind11::module module, pybind11::class_< SubspaceBasis, options... > &cls )
    {
      using RootBasis = typename SubspaceBasis::RootBasis;

      // Use default constructed TreePath. This requires, that
      // the PrefixPath template parameter of SubspaceBasis is
      // a fully static tree path.
      auto construct = [] ( const RootBasis &rootBasis ) { return new SubspaceBasis( rootBasis, {} ); };
      registerBasisType ( module, cls, construct, std::false_type{} );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_FUNCTIONS_SUBSPACEBASIS_HH
