// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localfiniteelementvariant.hh>
#include <dune/localfunctions/raviartthomas.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas0cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas0cube3d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas02d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas03d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube3d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas12d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas2cube2d.hh>

#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>

namespace Dune {
namespace Functions {

namespace Impl {

  template<int dim, typename D, typename R, std::size_t k>
  struct RaviartThomasSimplexLocalInfo
  {
    // Dummy type, must be something that we can have a std::unique_ptr to
    using FiniteElement = void*;
  };

  template<typename D, typename R>
  struct RaviartThomasSimplexLocalInfo<2,D,R,0>
  {
    using FiniteElement = RT02DLocalFiniteElement<D,R>;
  };

  template<typename D, typename R>
  struct RaviartThomasSimplexLocalInfo<2,D,R,1>
  {
    using FiniteElement = RT12DLocalFiniteElement<D,R>;
  };

  template<typename D, typename R>
  struct RaviartThomasSimplexLocalInfo<3,D,R,0>
  {
    using FiniteElement = RT03DLocalFiniteElement<D,R>;
  };

  template<int dim, typename D, typename R, std::size_t k>
  struct RaviartThomasCubeLocalInfo
  {
    // Dummy type, must be something that we can have a std::unique_ptr to
    using FiniteElement = void*;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<2,D,R,0>
  {
    using FiniteElement = RT0Cube2DLocalFiniteElement<D,R>;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<2,D,R,1>
  {
    using FiniteElement = RT1Cube2DLocalFiniteElement<D,R>;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<2,D,R,2>
  {
    using FiniteElement = RT2Cube2DLocalFiniteElement<D,R>;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<3,D,R,0>
  {
    using FiniteElement = RT0Cube3DLocalFiniteElement<D,R>;
  };

  template<typename D, typename R>
  struct RaviartThomasCubeLocalInfo<3,D,R,1>
  {
    using FiniteElement = RT1Cube3DLocalFiniteElement<D,R>;
  };

  template<int dim, typename D, typename R, std::size_t k>
  struct RaviartThomasPyramidLocalInfo
  {
    // Dummy type, must be something that we can have a std::unique_ptr to
    using FiniteElement = void*;
  };

  template<typename D, typename R>
  struct RaviartThomasPyramidLocalInfo<3,D,R,0>
  {
    using FiniteElement = RT0PyramidLocalFiniteElement<D,R>;
  };

  template<int dim, typename D, typename R, std::size_t k>
  struct RaviartThomasPrismLocalInfo
  {
    // Dummy type, must be something that we can have a std::unique_ptr to
    using FiniteElement = void*;
  };

  template<typename D, typename R>
  struct RaviartThomasPrismLocalInfo<3,D,R,0>
  {
    using FiniteElement = RT0PrismLocalFiniteElement<D,R>;
  };

  template<typename GV, int dim, typename R, std::size_t k>
  class RaviartThomasLocalFiniteElementMap
  {
    using D = typename GV::ctype;
    constexpr static bool hasFixedElementType = Capabilities::hasSingleGeometryType<typename GV::Grid>::v;

    using CubeFiniteElement    = typename RaviartThomasCubeLocalInfo<dim, D, R, k>::FiniteElement;
    using SimplexFiniteElement = typename RaviartThomasSimplexLocalInfo<dim, D, R, k>::FiniteElement;
    using PyramidFiniteElement = typename RaviartThomasPyramidLocalInfo<dim, D, R, k>::FiniteElement;
    using PrismFiniteElement = typename RaviartThomasPrismLocalInfo<dim, D, R, k>::FiniteElement;

  public:

    using T = LocalBasisTraits<D, dim, FieldVector<D,dim>, R, dim, FieldVector<R,dim>, FieldMatrix<D,dim,dim> >;

    constexpr static unsigned int  topologyId = Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId;  // meaningless if hasFixedElementType is false
    constexpr static GeometryType type = GeometryType(topologyId, GV::dimension);

    using FiniteElement = std::conditional_t<hasFixedElementType,
                                           std::conditional_t<type.isCube(), CubeFiniteElement, SimplexFiniteElement>,
                              std::conditional_t<dim == 3,
                                           LocalFiniteElementVariant<CubeFiniteElement, SimplexFiniteElement, PyramidFiniteElement, PrismFiniteElement>,
                                           LocalFiniteElementVariant<CubeFiniteElement, SimplexFiniteElement> > >;

    // Each element facet can have its orientation reversed, hence there are
    // 2^#facets different variants.
    static std::size_t numVariants(GeometryType type)
    {
      auto numFacets = referenceElement<D,dim>(type).size(1);
      return power(2,numFacets);
    }

    RaviartThomasLocalFiniteElementMap(const GV& gv)
      : elementMapper_(gv, mcmgElementLayout())
    {
      update(gv);
    }

    void update (const GV& gv)
    {
      elementMapper_.update(gv);
      if constexpr (hasFixedElementType)
      {
        variants_.resize(numVariants(type));
        for (size_t i = 0; i < numVariants(type); i++)
          variants_[i] = FiniteElement(i);
      }
      else
      {
        // for mixed grids add offset for cubes
        size_t numVariantsSimplex = numVariants(GeometryTypes::simplex(dim));
        size_t numVariantsCube = numVariants(GeometryTypes::cube(dim));

        variants_.resize(numVariantsSimplex + numVariantsCube);
        for (size_t i = 0; i < numVariantsSimplex; i++)
          variants_[i] = SimplexFiniteElement(i);
        for (size_t i = 0; i < numVariantsCube; i++)
          variants_[i + numVariantsSimplex] = CubeFiniteElement(i);
        if constexpr (dim == 3)
        {
          size_t numVariantsPyramid = numVariants(GeometryTypes::pyramid);
          size_t numVariantsPrism = numVariants(GeometryTypes::prism);

          variants_.resize(numVariantsSimplex + numVariantsCube + numVariantsPyramid + numVariantsPrism );

          for (size_t i = 0; i < numVariantsPyramid; i++)
            variants_[i + numVariantsSimplex + numVariantsCube] = PyramidFiniteElement(i);
          for (size_t i = 0; i < numVariantsPrism; i++)
            variants_[i + numVariantsSimplex + numVariantsCube + numVariantsPyramid] = PrismFiniteElement(i);
        }
      }

      orient_.resize(gv.size(0));
      for(const auto& cell : elements(gv))
      {
        unsigned int myId = elementMapper_.index(cell);
        orient_[myId] = 0;

        for (const auto& intersection : intersections(gv,cell))
        {
          if (intersection.neighbor() && (gv.contains(intersection.outside())) && (elementMapper_.index(intersection.outside()) > myId))
            orient_[myId] |= (1 << intersection.indexInInside());
        }

        // for mixed grids add offset for cubes and pyramids
        if constexpr (!hasFixedElementType)
        {
          size_t numVariantsSimplex = numVariants(GeometryTypes::simplex(dim));
          size_t numVariantsCube = numVariants(GeometryTypes::cube(dim));

          if (cell.type().isCube())
            orient_[myId] += numVariantsSimplex;
          if constexpr (dim == 3)
          {
            size_t numVariantsPyramid = numVariants(GeometryTypes::pyramid);

            if (cell.type().isPyramid())
              orient_[myId] += numVariantsSimplex + numVariantsCube;
            if (cell.type().isPrism())
              orient_[myId] += numVariantsSimplex + numVariantsCube + numVariantsPyramid;
          }
        }
      }
    }

    template<class EntityType>
    const FiniteElement& find(const EntityType& e) const
    {
      return variants_[orient_[elementMapper_.index(e)]];
    }

    private:
      std::vector<FiniteElement> variants_;
      Dune::MultipleCodimMultipleGeomTypeMapper<GV> elementMapper_;
      std::vector<unsigned char> orient_;
  };


} // namespace Impl


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   RaviartThomasPreBasis
//   RaviartThomasNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k>
class RaviartThomasNode;

template<typename GV, int k>
class RaviartThomasPreBasis :
  public LeafPreBasisMapperMixin<GV>
{
  using Base = LeafPreBasisMapperMixin<GV>;

  // the layout is defined in terms of a MCMGLayout specialized for k == 0 or 1
  static MCMGLayout dofLayout()
  {
    return [](GeometryType gt, int gridDim) -> size_t {
      if ((gt.isPyramid()) and (k==0))
        return 1;
      if (gt.dim() == gridDim)
        return gt.isCube() ? ((dim == 2) ? k*(k+1)*dim : k*(k+1)*(k+1)*dim) : k*dim;
      if (gt.dim() == gridDim-1)
        return gt.isCube() ? (dim-2)*2*k+k+1 : (dim-1)*k+1 ;
      return 0;
    };
  }

  static const int dim = GV::dimension;
  using FiniteElementMap = typename Impl::RaviartThomasLocalFiniteElementMap<GV, dim, double, k>;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  using Node = RaviartThomasNode<GV, k>;

  /** \brief Constructor for a given grid view object */
  RaviartThomasPreBasis(const GridView& gv) :
    Base(gv, dofLayout()),
    finiteElementMap_(gv)
  {
    // Currently there are some unresolved bugs with hybrid grids and higher order Raviart-Thomas elements
    if (gv.indexSet().types(0).size() > 1 and k>0)
      DUNE_THROW(Dune::NotImplemented, "Raviart-Thomas basis with index k>0 is only implemented for grids with a single element type");

    for(auto type : gv.indexSet().types(0))
      if (!type.isSimplex() && !type.isCube() && !type.isPyramid() && !type.isPrism())
        DUNE_THROW(Dune::NotImplemented, "Raviart-Thomas elements are only implemented for grids with simplex, cube, pyramid or prism elements.");
  }

  /** \brief Update the stored grid view, to be called if the grid has changed */
  void update(const GridView& gv)
  {
    Base::update(gv);
    finiteElementMap_.update(gv);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{&finiteElementMap_};
  }

protected:
  FiniteElementMap finiteElementMap_;
};



template<typename GV, int k>
class RaviartThomasNode :
  public LeafBasisNode
{
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElementMap = typename Impl::RaviartThomasLocalFiniteElementMap<GV, dim, double, k>;
  using FiniteElement = Impl::GlobalValuedLocalFiniteElement<Impl::ContravariantPiolaTransformator,
                                                             typename FiniteElementMap::FiniteElement,
                                                             Element>;

  RaviartThomasNode(const FiniteElementMap* finiteElementMap) :
    element_(nullptr),
    finiteElementMap_(finiteElementMap)
  { }

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_.bind((finiteElementMap_->find(*element_)), e);
    this->setSize(finiteElement_.size());
  }

protected:

  FiniteElement finiteElement_;
  const Element* element_;
  const FiniteElementMap* finiteElementMap_;
};

namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a Raviart-Thomas pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k Order of the Raviart-Thomas element
 */
template<std::size_t k>
auto raviartThomas()
{
  return [](const auto& gridView) {
    return RaviartThomasPreBasis<std::decay_t<decltype(gridView)>, k>(gridView);
  };
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a k-th-order Raviart Thomas finite element space
 *
 * TODO: Fix this for grids with more than one element type
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using RaviartThomasBasis = DefaultGlobalBasis<RaviartThomasPreBasis<GV, k> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASBASIS_HH
