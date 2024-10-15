// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_MAPPERUTILITIES_HH
#define DUNE_FUNCTIONS_COMMON_MAPPERUTILITIES_HH

#include <vector>

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



} // end namespace Impl

} // end namespace Dune::Functions


#endif // end namespace DUNE_FUNCTIONS_COMMON_MAPPERUTILITIES_HH
