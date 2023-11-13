// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICALLAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICALLAGRANGEBASIS_HH

#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/hierarchical/hierarchicalp2.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lfeprebasismixin.hh>

#include <dune/geometry/type.hh>

namespace Dune {
  namespace Functions {

    // *****************************************************************************
    // Implementation for Hierarchical Lagrange Basis
    //
    // - only orders k=1,2 are implemented up to now
    // - order k=1 is identical to the standard Lagrange basis
    // - implementation is restricted to simplex grids
    //
    // *****************************************************************************

    /**
     * \brief A pre-basis for a hierarchical Lagrange basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV  The grid view that the FE basis is defined on
     * \tparam k   The polynomial order of ansatz functions, (0 < k <= 2)
     * \tparam R   Range field-type used for shape function values
     */
    template<typename GV, int k, typename R = double>
    class HierarchicalLagrangePreBasis;

    template<typename GV, typename R>
    class HierarchicalLagrangePreBasis<GV,1,R>
        : public LFEPreBasisMixin<GV, LagrangeSimplexLocalFiniteElement<typename GV::ctype,R,GV::dimension,1>>
    {
      using Base = LFEPreBasisMixin<GV, LagrangeSimplexLocalFiniteElement<typename GV::ctype,R,GV::dimension,1>>;
    public:
      HierarchicalLagrangePreBasis (const GV& gridView) :
        Base(gridView, [](GeometryType gt, int) -> std::size_t { return gt.isVertex() ? 1 : 0; })
      {
        for (auto gt : gridView.indexSet().types(0)) {
          if (!gt.isSimplex())
            DUNE_THROW(Dune::NotImplemented,
              "Hierarchical Lagrange basis only implemented for simplex grids.");
        }
      }
    };

    template<typename GV, typename R>
    class HierarchicalLagrangePreBasis<GV,2,R>
        : public LFEPreBasisMixin<GV, HierarchicalP2LocalFiniteElement<typename GV::ctype,R,GV::dimension>>
    {
      using Base = LFEPreBasisMixin<GV, HierarchicalP2LocalFiniteElement<typename GV::ctype,R,GV::dimension>>;
    public:
      HierarchicalLagrangePreBasis (const GV& gridView) :
        Base(gridView, [](GeometryType gt, int) -> std::size_t { return (gt.isVertex() || gt.isLine()) ? 1 : 0; })
      {
        for (auto gt : gridView.indexSet().types(0)) {
          if (!gt.isSimplex())
            DUNE_THROW(Dune::NotImplemented,
              "Hierarchical Lagrange basis only implemented for simplex grids.");
        }
      }
    };


    namespace BasisFactory {

      /**
       * \brief A factory that can create a HierarchicalLagrange pre-basis
       *
       * \ingroup FunctionSpaceBasesImplementations
       *
       * \tparam k   The polynomial order of the ansatz functions (0 < k <= 2)
       * \tparam R   The range field-type of the local basis
       */
      template<int k, typename R=double>
      auto hierarchicalLagrange()
      {
        static_assert(0 < k && k <= 2);
        return [](const auto& gridView) {
          return HierarchicalLagrangePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
        };
      }

    } // end namespace BasisFactory

    /** \brief Basis of a scalar Hierarchical Lagrange finite element space
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV The GridView that the space is defined on
     * \tparam k The order of the basis (0 < k <= 2)
     * \tparam R The range type of the local basis
     *
     *  \note currently only supports simplex grids
     */
    template<typename GV, int k, typename R=double>
    using HierarchicalLagrangeBasis = DefaultGlobalBasis<HierarchicalLagrangePreBasis<GV, k, R> >;

  } // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICALLAGRANGEBASIS_HH
