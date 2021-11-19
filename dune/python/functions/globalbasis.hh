#ifndef DUNE_PYTHON_FUNCTIONS_GLOBALBASIS_HH
#define DUNE_PYTHON_FUNCTIONS_GLOBALBASIS_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

#include <dune/python/common/dimrange.hh>
#include <dune/python/common/fmatrix.hh>
#include <dune/python/common/fvector.hh>
#include <dune/python/functions/discretefunction.hh>
#include <dune/python/functions/interpolate.hh>
#include <dune/python/functions/tree.hh>

#include <dune/python/pybind11/complex.h>
#include <dune/python/pybind11/pybind11.h>


namespace Dune
{

  namespace Python
  {

    template <class Basis>
    struct LocalViewWrapper : public Basis::LocalView
    {
      typedef typename Basis::LocalView Base;
      typedef typename Base::Element EntityType;
      LocalViewWrapper(const Basis &b) : Base(b) {}

      std::vector<int> index(int idx) const
      {
        // call index in the base class
        auto ind = Base::index(idx);
        std::vector<int> ret(ind.size());
        for (int i=0;i<ind.size();++i) ret[i] = ind[i];
        return ret;
      }

      void bind ( pybind11::object &obj )
      {
        obj_ = obj;
        const EntityType &entity = obj.template cast<const EntityType&>();
        Base::bind(entity);
      }
      void unbind ( )
      {
        Base::unbind();
        obj_.release();
      }
      pybind11::object obj_;
    };

    template< class GlobalBasis, class... options >
    DUNE_EXPORT void registerGlobalBasis ( pybind11::module module, pybind11::class_< GlobalBasis, options... > &cls )
    {
      using pybind11::operator""_a;

      typedef Dune::TypeTree::HybridTreePath<> DefaultTreePath;

      const std::size_t dimRange = DimRange< typename GlobalBasis::PreBasis::Node >::value;

      cls.def( pybind11::init( [] ( const typename GlobalBasis::GridView &gridView ) { return new GlobalBasis( gridView ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def( "__len__", [](const GlobalBasis& self) { return self.dimension(); } );

      cls.def_property_readonly( "dimRange", [] ( pybind11::handle self ) { return pybind11::int_( dimRange ); } );
      cls.def_property( "gridView",
                        [](const GlobalBasis& basis) { return basis.gridView(); },
                        [](GlobalBasis& basis, const typename GlobalBasis::GridView& gridView) { basis.update(gridView); });

      typedef LocalViewWrapper< GlobalBasis > LocalView;
      auto lv = insertClass< LocalView >( module, "LocalView",
          GenerateTypeName("Dune::Python::LocalViewWrapper", MetaType<GlobalBasis>()),
          IncludeFiles{"dune/python/functions/globalbasis.hh"}).first;
      lv.def( "bind", &LocalView::bind );
      lv.def( "unbind", &LocalView::unbind );
      lv.def( "index", [] ( const LocalView &localView, int index ) { return localView.index( index ); });
      lv.def( "__len__", [] ( LocalView &self ) -> int { return self.size(); } );

      Functions::registerTree<typename LocalView::Tree>(lv);
      lv.def("tree", [](const LocalView& view) { return view.tree(); });

      cls.def( "localView", [] ( const GlobalBasis &self ) -> LocalView { return LocalView( self ); }, pybind11::keep_alive< 0, 1 >() );
      cls.def_property_readonly( "dimension", [] ( const GlobalBasis &self ) -> int { return self.dimension(); } );

      cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, double> );
      cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, bool> );
      cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, int> );

      typedef Dune::FieldVector< double, dimRange > Range;
      typedef Dune::Functions::DiscreteGlobalBasisFunction< GlobalBasis, HierarchicPythonVector< double >, DefaultNodeToRangeMap< GlobalBasis, DefaultTreePath >, Range > DiscreteFunction;
      auto clsDiscreteFunction = insertClass< DiscreteFunction >( module, "DiscreteFunction", GenerateTypeName( cls, "DiscreteFunction" ) );
      registerDiscreteFunction( module, clsDiscreteFunction.first );

      cls.def("asFunction", [] ( GlobalBasis &self, pybind11::buffer dofVector ) {
          auto nodeToRangeMapPtr =
            std::make_shared< const DefaultNodeToRangeMap<
              GlobalBasis, DefaultTreePath >
                              >(
                                makeDefaultNodeToRangeMap(self, DefaultTreePath()));
          std::shared_ptr<const GlobalBasis> basisPtr = Dune::wrap_or_move( self );
          auto vectorPtr = std::make_shared< const HierarchicPythonVector< double > >( dofVector );
          return new DiscreteFunction( basisPtr,
                                       vectorPtr,
                                       nodeToRangeMapPtr);
        }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >(), "dofVector"_a );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_FUNCTIONS_GLOBALBASIS_HH