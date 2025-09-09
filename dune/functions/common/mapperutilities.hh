// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_MAPPERUTILITIES_HH
#define DUNE_FUNCTIONS_COMMON_MAPPERUTILITIES_HH

#include <vector>
#include <bitset>

#include <dune/common/rangeutilities.hh>

#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune::Functions {

// All utilities are in Impl:: and thus considered implementation
// details for now. However, one may want to think about making
// them public. Then they could also be put into dune-grid,
// since there's nothing dune-function specific about them.
namespace Impl {



  // Helper class providing an unordered range
  // of global indices associated to the element
  // within a MultipleCodimMultipleGeomTypeMapper.
  // This has to be bound to an element, before
  // it can be used.
  template<class GV>
  class MapperElementSubIndices
  {
    using IndexContainer = std::vector<typename Dune::MultipleCodimMultipleGeomTypeMapper<GV>::Index>;
  public:
    using GridView = GV;
    using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using Index = typename Mapper::Index;
    using Element = typename GridView::template Codim<0>::Entity;

    MapperElementSubIndices(const Mapper& mapper)
      : mapper_(mapper)
    {}

    // Bind to given element and precompute all indices.
    void bind(const Element& element)
    {
      constexpr auto dimension = GridView::dimension;
      indices_.clear();
      auto referenceElement = Dune::referenceElement<double, dimension>(element.type());
      for (auto codim : Dune::range(dimension + 1)) {
        for (auto subEntity : Dune::range(referenceElement.size(codim))) {
          std::size_t c = mapper_.layout()(referenceElement.type(subEntity, codim), dimension);
          if (c > 0) {
            std::size_t firstIndex = mapper_.subIndex(element, subEntity, codim);
            for (auto j : Dune::range(firstIndex, firstIndex + c)) {
              indices_.push_back(j);
            }
          }
        }
      }
    }

    auto begin() const
    {
      return indices_.begin();
    }

    auto end() const
    {
      return indices_.end();
    }

  private:
    const Mapper mapper_;
    IndexContainer indices_;
  };



  // Helper function computing an average mesh size per subentity
  // by averaging over the adjacent elements. This only considers
  // the subentities handled by the given mapper and returns a
  // std::vector<FieldType> of mesh sizes indexed according to the mapper.
  //
  // The average is determined by first computing the average volume
  // of adjacent elements and then taking the d-th root for a d-dimensional
  // grid.
  //
  // This operation has O(n) runtime (with n=mapper.size()),
  // allocates O(n) memory for the returned vector and additional
  // O(n) temporary memory.
  template<class FieldType = double, class Mapper>
  auto computeAverageSubEntityMeshSize(const Mapper& mapper)
  {
    constexpr auto dimension = Mapper::GridView::dimension;
    std::vector<unsigned int> adjacentElements(mapper.size(), 0);
    std::vector<FieldType> subEntityMeshSize(mapper.size(), 0.0);
    auto subIndices = Impl::MapperElementSubIndices(mapper);
    for(const auto& element : Dune::elements(mapper.gridView()))
    {
      auto A = element.geometry().volume();
      subIndices.bind(element);
      for(auto i : subIndices)
      {
        subEntityMeshSize[i] += A;
        ++(adjacentElements[i]);
      }
    }
    for(auto i : Dune::range(mapper.size()))
      subEntityMeshSize[i] = std::pow(subEntityMeshSize[i]/adjacentElements[i], 1./dimension);
    return subEntityMeshSize;
  }

  /**
   * \brief Computes a consistent orientation of edges. The orientation is chosen,
   * such that the tangent of an edge points from the vertex with the lower global Id
   * to the vertex with the larger global Id. This can agree with the reference orientation or not.
   * \tparam ElementMapper  Mappertype providing the gridView and an index method
   * \param  mapper         The mapper to compute the storage order of orientations
   * \returns A vector with a std::bitset<3> for each element in the grid, order according to mapper.
   * The bitset contains one bit for each edge, which is O iff the global orientation agrees with the
   * reference orientation.
   * \note This only works for a two dimensional, purely simplicial, grid.
   */
  template<class ElementMapper>
  std::vector<std::bitset<3>> computeEdgeOrientations(ElementMapper mapper)
  {
    constexpr int dim = 2;
    static_assert(dim == ElementMapper::GridView::dimension);

    auto const& gridView = mapper.gridView();
    std::vector<std::bitset<3>> orientations;
    orientations.resize(gridView.size(0));
    // compute orientation for all elements
    auto const& idSet = gridView.grid().globalIdSet();

    for (const auto &element : elements(gridView))
    {
      const auto &refElement = referenceElement(element);
      auto elementIndex = mapper.index(element);

      std::bitset<3>& orientation = orientations[elementIndex];


      for (std::size_t i = 0; i < element.subEntities(dim - 1); i++)
      {
        // Local vertex indices within the element are ordered, localV0 < localV1
        auto localV0 = refElement.subEntity(i, dim - 1, 0, dim);
        auto localV1 = refElement.subEntity(i, dim - 1, 1, dim);

        // Global vertex indices within the grid
        auto globalV0 = idSet.subId(element, localV0, dim);
        auto globalV1 = idSet.subId(element, localV1, dim);

        // The edge is flipped if the local ordering disagrees with global ordering
        orientation[i] = globalV0 > globalV1;
      }
    }
    return orientations;
  }

} // end namespace Impl

} // end namespace Dune::Functions


#endif // end namespace DUNE_FUNCTIONS_COMMON_MAPPERUTILITIES_HH
