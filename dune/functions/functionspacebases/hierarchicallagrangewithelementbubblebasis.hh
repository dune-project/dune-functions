// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HierarchicalLagrangeWOTHELEMENTBUBBLEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HierarchicalLagrangeWOTHELEMENTBUBBLEBASIS_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

#include <dune/localfunctions/hierarchical/hierarchicalp1withelementbubble.hh>
#include <dune/localfunctions/hierarchical/hierarchicalp2withelementbubble.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lfeprebasismixin.hh>

#include <dune/geometry/type.hh>

namespace Dune {
  namespace Functions {

    // *****************************************************************************
    // Implementation for Hierarchical Lagrange Basis (of order 1-2) with an
    // element bubble function (of order dim+1).
    //
    // - currently only supports simplex grids
    //
    // *****************************************************************************

    /**
     * \brief A pre-basis for a hierarchical basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV  The grid view that the FE basis is defined on
     * \tparam k   The polynomial order of ansatz functions (0 < k <= 2)
     * \tparam R   Range field-type used for shape function values
     */
    template<typename GV, int k, typename R=double>
    class HierarchicalLagrangeWithElementBubblePreBasis;

    template<typename GV, typename R>
    class HierarchicalLagrangeWithElementBubblePreBasis<GV,1,R>
        : public LFEPreBasisMixin<GV, HierarchicalP1WithElementBubbleLocalFiniteElement<typename GV::ctype,R,GV::dimension>>
    {
      using Base = LFEPreBasisMixin<GV, HierarchicalP1WithElementBubbleLocalFiniteElement<typename GV::ctype,R,GV::dimension>>;
    public:
      HierarchicalLagrangeWithElementBubblePreBasis (const GV& gridView) :
        Base(gridView, [](GeometryType gt, int) { return (gt.isVertex() || gt.dim() == GV::dimension) ? 1 : 0; })
      {
        for (auto gt : gridView.indexSet().types(0)) {
          if (!gt.isSimplex())
            DUNE_THROW(Dune::NotImplemented,
              "Hierarchical Lagrange basis only implemented for simplex grids.");
        }
      }
    };

    template<typename GV, typename R>
    class HierarchicalLagrangeWithElementBubblePreBasis<GV,2,R>
        : public LFEPreBasisMixin<GV, HierarchicalP2WithElementBubbleLocalFiniteElement<typename GV::ctype,R,GV::dimension>>
    {
      using Base = LFEPreBasisMixin<GV, HierarchicalP2WithElementBubbleLocalFiniteElement<typename GV::ctype,R,GV::dimension>>;
    public:
      HierarchicalLagrangeWithElementBubblePreBasis (const GV& gridView) :
        Base(gridView, [](GeometryType gt, int) { return (gt.dim() <= 1 || gt.dim() == GV::dimension) ? 1 : 0; })
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
       * \brief A factory that can create a HierarchicalLagrangeWithElementBubble pre-basis
       *
       * \ingroup FunctionSpaceBasesImplementations
       *
       * \tparam k   The polynomial order of the ansatz functions (0 < k <= 2)
       * \tparam R   The range field-type of the local basis
       */
      template<int k, typename R=double>
      auto hierarchicalLagrangeWithElementBubble()
      {
        static_assert(0 < k && k <= 2);
        return [](const auto& gridView) {
          return HierarchicalLagrangeWithElementBubblePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
        };
      }

      //! Explicit factory for k=1 for the HierarchicalLagrangeWithElementBubblePreBasis pre-basis
      template<typename R=double>
      auto hierarchicalP1B()
      {
        return hierarchicalLagrangeWithElementBubble<1,R>();
      }

      //! Explicit factory for k=2 for the HierarchicalLagrangeWithElementBubblePreBasis pre-basis
      template<typename R=double>
      auto hierarchicalP2B()
      {
        return hierarchicalLagrangeWithElementBubble<2,R>();
      }

    } // end namespace BasisFactory

    /** \brief Basis of a Hierarchical Lagrange finite element space with element bubble functions
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV The GridView that the space is defined on
     * \tparam k The order of the basis (0 < k <= 2)
     * \tparam R The range field-type of the local basis
     *
     *  \note currently only supports simplex grids
     */
    template<typename GV, int k, typename R=double>
    using HierarchicalLagrangeWithElementBubbleBasis = DefaultGlobalBasis<HierarchicalLagrangeWithElementBubblePreBasis<GV, k, R> >;

  } // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HierarchicalLagrangeWOTHELEMENTBUBBLEBASIS_HH
