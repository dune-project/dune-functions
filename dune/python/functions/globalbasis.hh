#ifndef DUNE_PYTHON_FUNCTIONS_GLOBALBASIS_HH
#define DUNE_PYTHON_FUNCTIONS_GLOBALBASIS_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/classname.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

#include <dune/typetree/nodetags.hh>

#include <dune/python/common/dimrange.hh>
#include <dune/python/common/fmatrix.hh>
#include <dune/python/common/fvector.hh>
#include <dune/python/functions/discretefunction.hh>
#include <dune/python/functions/interpolate.hh>
#include <dune/python/functions/tree.hh>

#include <dune/python/pybind11/complex.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace PYBIND11_NAMESPACE
{
  namespace detail
  {
    /** \brief Implement type_caster for Dune::ReservedVector
     *
     * This uses the standard array_caster implementation for resizeable types.
     *
     * \warning array_caster is not documented officially. But as long as we have
     * a copy of pybind11 vendored in dune-common using it should be fine.
     */
    template<class T, int n>
    struct type_caster<Dune::ReservedVector<T, n> > : array_caster<Dune::ReservedVector<T, n>, T, true> {};
  }
}

namespace Dune
{

  namespace Python
  {

    namespace detail {

      // specialization of the DimRange utility from dune-common
      template <class T>
      struct DimRange<T, std::enable_if_t<std::is_same_v<typename T::NodeTag, Dune::TypeTree::CompositeNodeTag>> >
        : public DimRange<typename T::ChildTypes>
      {};

      template <class T>
      struct DimRange<T, std::enable_if_t<std::is_same_v<typename T::NodeTag, Dune::TypeTree::PowerNodeTag>> >
        : public std::integral_constant<std::size_t, T::degree() * DimRange<typename T::ChildType>::value>
      {};

      template <class T>
      struct DimRange<T, std::enable_if_t<std::is_same_v<typename T::NodeTag, Dune::TypeTree::LeafNodeTag>> >
        : public std::integral_constant<std::size_t, T::FiniteElement::Traits::LocalBasisType::Traits::dimRange>
      {};

    } // end namespace detail


    template <class Basis>
    struct LocalViewWrapper : public Basis::LocalView
    {
      typedef typename Basis::LocalView Base;
      typedef typename Base::Element EntityType;
      LocalViewWrapper(const Basis &b) : Base(b) {}
      LocalViewWrapper(const Base &lv) : Base(lv) {}

      Base& base()
      {
        return *this;
      }

      const Base& base() const
      {
        return *this;
      }

      std::vector<int> index(int idx) const
      {
        // call index in the base class
        auto ind = base().index(idx);
        std::vector<int> ret(ind.size());
        for (std::size_t i=0;i<ind.size();++i) ret[i] = ind[i];
        return ret;
      }

      void bind ( pybind11::object &obj )
      {
        obj_ = obj;
        const EntityType &entity = obj.template cast<const EntityType&>();
        base().bind(entity);
      }
      void unbind ( )
      {
        base().unbind();
        obj_.release();
      }
      pybind11::object obj_;
    };

    // The type of the value of a global finite element function at a particular point
    template<typename K, unsigned int n>
    struct RangeType
    {
      using type = Dune::FieldVector< K, n >;
      static void registerRange(pybind11::module scope)
      {
        registerFieldVector<K,n>(scope);
      }
    };

    template<typename K>
    struct RangeType<K,1>
    {
      using type = K;
      static void registerRange(pybind11::module scope) {} // nothing to register, as K is a basic type
    };

    template< class GlobalBasis, class... options, class ConstructCall, bool hasUpdate>
    void registerBasisType ( pybind11::module module, pybind11::class_< GlobalBasis, options... > &cls, ConstructCall constructCall, std::bool_constant<hasUpdate>)
    {
      using pybind11::operator""_a;
      using GridView = typename GlobalBasis::GridView;

      const std::size_t dimRange = DimRange< typename GlobalBasis::LocalView::Tree >::value;
      const std::size_t dimWorld = GridView::dimensionworld;

      cls.def( pybind11::init( constructCall ), pybind11::keep_alive< 1, 2 >() );
      cls.def( "__len__", [](const GlobalBasis& self) { return self.dimension(); } );
      cls.def("size",
               pybind11::overload_cast<>(&GlobalBasis::size, pybind11::const_),
               "Return number of possible values for next position in empty multi index" );
      cls.def("size",
              pybind11::overload_cast<const typename GlobalBasis::SizePrefix&>(&GlobalBasis::size, pybind11::const_),
              "Return number of possible values for next position in multi index");

      cls.def_property_readonly( "dimRange", [] ( pybind11::handle self ) { return pybind11::int_( dimRange ); } );

      if constexpr (hasUpdate)
      {
        cls.def_property( "gridView",
                          [](const GlobalBasis& basis) { return basis.gridView(); },
                          [](GlobalBasis& basis, const GridView& gridView) { basis.update(gridView); });
      }
      else
      {
        cls.def_property_readonly( "gridView", [](const GlobalBasis& basis) { return basis.gridView(); });
      }

      using LocalView = LocalViewWrapper< GlobalBasis >;
      using Tree = typename LocalView::Tree;

      auto includes = IncludeFiles{"dune/python/functions/globalbasis.hh"};
      auto lv = insertClass< LocalView >( module, "LocalView",
          GenerateTypeName("Dune::Python::LocalViewWrapper", MetaType<GlobalBasis>()),
          includes).first;
      lv.def( "bind", &LocalView::bind );
      lv.def( "unbind", &LocalView::unbind );
      lv.def( "element", &LocalView::element );
      lv.def( "index", [] ( const LocalView &localView, int index ) { return localView.index( index ); });
      lv.def( "__len__", [] ( LocalView &self ) -> int { return self.size(); } );
      lv.def( "size", [] ( LocalView &self ) -> int { return self.size(); } );
      lv.def( "maxSize", [] ( LocalView &self ) -> int { return self.maxSize(); } );

      Functions::registerTree<Tree>(lv);
      lv.def("tree", [](const LocalView& view) { return Dune::stackobject_to_shared_ptr(view.tree()); });

      cls.def( "localView", [] ( const GlobalBasis &self ) -> LocalView { return self.localView(); }, pybind11::keep_alive< 0, 1 >() );
      cls.def_property_readonly( "dimension", [] ( const GlobalBasis &self ) -> int { return self.dimension(); } );

      // Register the 'interpolate' method
      // TODO: Currently we only register the method for cases where either scalars or FieldVector<double,dimRange>
      // are reasonable value types for functions to be interpolated. This excludes various not-very-exotic
      // composite bases.  Support for them is planned for the future.
      if constexpr (dimRange==1)
      {
        cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, double> );
        cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, bool> );
        cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, int> );
      }
      else if constexpr (GlobalBasis::LocalView::Tree::isLeaf or GlobalBasis::LocalView::Tree::isPower)
      {
        cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, double, FieldVector<double,dimRange> > );
        cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, bool, std::array<bool,dimRange> > );
        cls.def( "interpolate", &Dune::Python::Functions::interpolate<GlobalBasis, int, FieldVector<int,dimRange> > );
      }

      // Register various grid function types
      // TODO: As we do not currently have support for nested range types,
      // we register grid function types only for 'simple' bases.
      // A more general implementation is planned for the future.
      if constexpr (GlobalBasis::LocalView::Tree::isLeaf or GlobalBasis::LocalView::Tree::isPower)
      {
        using Range = typename RangeType< double, dimRange >::type;
        RangeType< double, dimRange >::registerRange(module);
        using Domain = Dune::FieldVector< double, dimWorld >;
        registerFieldVector<double,dimWorld>(module);
        using DiscreteFunction = Dune::Functions::DiscreteGlobalBasisFunction< GlobalBasis, HierarchicPythonVector< double >, Dune::Functions::HierarchicNodeToRangeMap, Range >;
        // register the HierarchicPythonVector
        Dune::Python::addToTypeRegistry<HierarchicPythonVector<double>>(
          GenerateTypeName("Dune::Python::HierarchicPythonVector", MetaType<double>()),
          {"dune/python/functions/discretefunction.hh"}
          );
        // and add the DiscreteFunction to our module
        auto clsDiscreteFunction = insertClass< DiscreteFunction >( module, "DiscreteFunction",
          GenerateTypeName( "Dune::Functions::DiscreteGlobalBasisFunction",
            MetaType<GlobalBasis>(),
            MetaType<HierarchicPythonVector< double >>(),
            "Dune::Functions::HierarchicNodeToRangeMap",
            MetaType<Range>()
            ), includes);
        // register the GridViewFunction and register the implicit conversion
        Dune::Python::addToTypeRegistry<Range(Domain)>(GenerateTypeName(className<Range(Domain)>()));
        using GridViewFunction = Dune::Functions::GridViewFunction<Range(Domain), GridView>;
        auto clsGridViewFunction = insertClass< GridViewFunction >( module, "GridViewFunction",
          GenerateTypeName( "Dune::Functions::GridViewFunction",
            MetaType<Range(Domain)>(),
            MetaType<GridView>()
            ), includes);
        clsGridViewFunction.first.def(pybind11::init<DiscreteFunction>());
        pybind11::implicitly_convertible<DiscreteFunction, GridViewFunction>();

        registerDiscreteFunction<GlobalBasis>( module, clsDiscreteFunction.first );

        cls.def("asFunction", [] ( GlobalBasis &self, pybind11::buffer dofVector ) {
            return new DiscreteFunction( self, HierarchicPythonVector<double>(dofVector), Dune::Functions::HierarchicNodeToRangeMap());
          }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >(), "dofVector"_a );
      }
    }

    template< class GlobalBasis, class... options >
    DUNE_EXPORT void registerGlobalBasis ( pybind11::module module, pybind11::class_< GlobalBasis, options... > &cls )
    {
      using GridView = typename GlobalBasis::GridView;
      auto construct = [] ( const GridView &gridView ) { return new GlobalBasis( gridView ); };
      registerBasisType ( module, cls, construct, std::true_type{} );
    }


  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_FUNCTIONS_GLOBALBASIS_HH
