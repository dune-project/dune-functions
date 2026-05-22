// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/typelist.hh>
#include <dune/common/indices.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismixin.hh>

#include <dune/grid/common/capabilities.hh>

namespace Dune {
namespace Functions {



namespace Experimental {

/**
 * \brief A class encoding the orientation of the subentities of an entity
 *
 * The FaceOrientations class provides the orientation of all faces
 * of a grid element in 2d and 3d relative to their canonical orientation
 * provided by the respective ReferenceElement.
 *
 * Since 0d-faces (a.k.a. vertices) only have a single orientation this is only
 * relevant for 1d-faces (edges) and 2d faces (triangles and quadrilaterals).
 * Notice that the orientations of a face are the automorphisms of the respective
 * edge graph. I.e. the orientations of an edge are elements from the
 * [symmetric group](https://en.wikipedia.org/wiki/Symmetric_group)
 * \f$S_2\f$ while the orientations of triangles and quadrilaterals are
 * elements from the [dihedral groups](https://en.wikipedia.org/wiki/Dihedral_group)
 * \f$D_3\f$ and \f$D_4\f$, respectively.
 * The internal encoding of a face orientation describes how its _global orientation_
 * can be mapped to its _local orientation_. These two orientations are defined as follows:
 *
 * * The `ReferenceElement` of a grid element provides a _canonical local orientation_
 *   for each face which is induced by the order of the face vertices in the `ReferenceElement`.
 *   This orientation associates zero-based consecutive local indices to the vertices of
 *   a face. Since the local orientation depends on the mapping of the element to the
 *   `ReferenceElement`, the same face may have different local orientations when viewed
 *   from different adjacent grid elements.
 *
 * * Given a set of globally unique vertex indices or ids, we can also associate a
 *   _unique global orientation_ to each face as follows: We first select the face
 *   vertex with the smallest global id and denote it by \f$P_0\f$. Among its neighboring
 *   vertices within the face (one for an edge, two for a 2d-face) we again select
 *   the one with the smallest global id and denote it by \f$P_1\f$. Then there is a unique
 *   automorphisms of the face that maps the edge \f$(P_0,P_1)\f$ to the edge with local
 *   vertex indices 0 and 1. We identify this automorphism with the face orientation.
 *
 * **Binary encoding of all face orientations**
 *
 * The FaceOrientations class now encodes these automorphisms for all
 * faces of the element in binary form. Since we have \f$|S_2|=2\f$ we can
 * store the edge orientation of each edge in a single bit, while we need
 * three bits to identify one of the \f$|D_3|=6\f$ or \f$|D_4|=8\f$ orientations
 * of a 2d-face. To furthermore distinguish between triangle orientations
 * from \f$D_3\f$ and quadrilateral orientations from \f$D_4\f$ we use another
 * bit indicating if the respective face is a triangle or quadrilateral.
 * Since we have at most 12 edges and 6 2d-faces for elements up to dimension 3
 * we can use a bitfield with 36 bits as underlying storage.
 * The 12 lowest order bits encode the 1d-face orientations. These are followed
 * by 6 block with 4 bits for the 2d-faces.
 *
 * * For a 1d-face (edge) the bit is unset if local and global
 *   orientation coincide and set if the edge needs to be flipped.
 * * For a 2d-face (triangle or quadrilateral) the global orientation
 *   can be mapped to the local one by first rotating the vertices
 *   counter-clockwise and then reflecting the face across the
 *   xy-diagonal. The two lowest order bits store the number
 *   of counter clockwise rotations and the third is set if a subsequent
 *   reflection is required. This encoding makes use of the fact
 *   that a generator for the groups \f$D_3\f$ and \f$D_4\f$ is given by
 *   \f$\{r,s\}\f$ where \f$r\f$ is an elementary rotation and \f$s\f$ a reflection.
 *   In the binary representation of the automorphisms \f$a=s^j r^i\f$
 *   the first two bits encode the exponent \f$i\f$ and the third
 *   one the exponent \f$j\f$.
 *   The fourth bit is unset for a triangle and set for a quadrilateral.
 *   Interpreting the bits from right to left we get:
 *   \f\[
 *      (\underbrace{b_4}_{\text{is quadrilateral}}, \underbrace{b_3}_{j = b_3},\underbrace{b_2,b_1}_{i = 2 b_2 + b_1})
 *   \f\]
 *
 * The 1 or 4 bits stored for each 1d- or 2d-face can be interpreted as
 * an index from 0,1 or 0,...,15, respectively.
 * This index identifies the orientation of the corresponding face and
 * is unique when viewing the face from different elements. It can e.g.
 * be used for a table lookup of the face DOF permutation required to
 * uniquely identify face DOFs for higher order Lagrange elements.
 *
 * The following tables depict all possible face orientations for triangle or
 * quadrilateral 2d-faces. Vertices and edges of the face are both enumerated
 * according to the `ReferenceElement`.
 * The global vertex ids are for simplicity denoted as 0,1,2,3 while any selection
 * of pairwise different less-than comparable ids can be used.
 * The primary edge is the one that induces the global orientations, i.e., connecting
 * the vertex with smallest id to its neighbors with smallest id. The edge orientations
 * in the table indicate if the (globally oriented) edge of the face is flipped
 * with respect to the `ReferenceElement`.
 *
 * **Triangle orientations**
 *
 * | Global vertex ids              | edge orientations | primary edge | rotations | reflection  | automorphism  | is quadrilateral  | bit-encoding  | index |
 * |--------------------------------|-------------------|--------------|-----------|-------------|---------------|-------------------|---------------|-------|
 * | 0,1,2                          | 0,0,0             | 0->1         | 0         | no          | \f$s^0 r^0\f$ | no                | 0000          | 0     |
 * | 1,2,0                          | 0,1,1             | 2->0         | 1         | no          | \f$s^0 r^1\f$ | no                | 0001          | 1     |
 * | 2,0,1                          | 1,1,0             | 1->0         | 2         | no          | \f$s^0 r^2\f$ | no                | 0010          | 2     |
 * | unused                         |                   |              | 3         | no          | \f$s^0 r^3\f$ | no                | 0011          | 3     |
 * | 0,2,1                          | 0,0,1             | 0->2         | 0         | yes         | \f$s^1 r^0\f$ | no                | 0100          | 4     |
 * | 2,1,0                          | 1,1,1             | 2->1         | 1         | yes         | \f$s^1 r^1\f$ | no                | 0101          | 5     |
 * | 1,0,2                          | 1,0,0             | 1->0         | 2         | yes         | \f$s^1 r^2\f$ | no                | 0110          | 6     |
 * | unused                         |                   |              | 3         | yes         | \f$s^1 r^3\f$ | no                | 0111          | 7     |
 *
 * Notice that using 4 bits is redundant for a triangle since only 0,1, or 2 rotations
 * are relevant leading to the unused indices 3 and 7. Similarly only six combinations
 * of edge orientations are possible for triangle
 *
 * **Quadrilateral orientations**
 *
 * | Global vertex ids              | edge orientations | primary edge | rotations | reflection  | automorphism  | is quadrilateral  | bit-encoding  | index |
 * |--------------------------------|-------------------|--------------|-----------|-------------|---------------|-------------------|---------------|-------|
 * | 0,1,2,3                        | 0,0,0,0           | 0->1         | 0         | no          | \f$s^0 r^0\f$ | yes               | 1000          | 8     |
 * | 0,1,3,2                        | 0,0,0,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 0,2,3,1                        | 0,1,0,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 1,2,0,3                        | 1,0,0,0           | 2->0         | 1         | no          | \f$s^0 r^1\f$ | yes               | 1001          | 9     |
 * | 1,3,0,2                        | 1,1,0,0           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 2,1,0,3                        | 1,0,1,0           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 2,3,1,0                        | 1,1,0,1           | 3->2         | 2         | no          | \f$s^0 r^2\f$ | yes               | 1010          | 10    |
 * | 3,2,1,0                        | 1,1,1,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 1,3,2,0                        | 0,1,0,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 2,0,3,1                        | 0,0,1,1           | 1->3         | 3         | no          | \f$s^0 r^3\f$ | yes               | 1011          | 11    |
 * | 3,0,2,1                        | 1,0,1,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 3,0,1,2                        | 1,0,1,0           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 0,2,1,3                        | 0,0,0,0           | 0->2         | 0         | yes         | \f$s^1 r^0\f$ | yes               | 1100          | 12    |
 * | 0,3,1,2                        | 0,1,0,0           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 0,3,2,1                        | 0,1,0,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 2,3,0,1                        | 1,1,0,0           | 2->3         | 1         | yes         | \f$s^1 r^1\f$ | yes               | 1101          | 13    |
 * | 3,2,0,1                        | 1,1,1,0           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 3,1,0,2                        | 1,0,1,0           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 2,1,3,0                        | 0,1,1,1           | 3->1         | 2         | yes         | \f$s^1 r^2\f$ | yes               | 1110          | 14    |
 * | 3,1,2,0                        | 1,1,1,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 1,2,3,0                        | 0,1,0,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 1,0,2,3                        | 0,0,1,0           | 1->0         | 3         | yes         | \f$s^1 r^3\f$ | yes               | 1111          | 15    |
 * | 1,0,3,2                        | 0,0,1,1           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 * | 2,0,1,3                        | 1,0,1,0           | ^            | ^         | ^           | ^             | ^                 | ^             | ^     |
 *
 * Notice that there are always three possible permutations of vertex ids for a
 * quadrilateral that lead to the same primary edge and thus the same orientation.
 * This corresponds to the fact that the dihedral group \f$D_4\f$ of quadrilateral
 * automorphisms is a subgroup of the symmetric group \f$S_4\f$ of all vertex permutations.
 * Furthermore for a quadrilateral the face orientation
 * cannot be deduced from the edge orientations alone.
 *
 * **Complexity bounds**
 *
 * The algorithm to compute the orientations of the faces uses in the constructor
 * of \ref FaceOrientations tries to minimize the number of vertex index comparisons
 * while preventing the use of full lookup tables
 * because the latter can become very large in 3d. In general we need one comparison
 * per edge and additional one or two comparisons per quadrilateral face.
 * The following table lists the number of minimal an maximal number of comparisons
 * of the implemented algorithm for the respective geometry types.
 * For comparison this also lists the number of comparisons
 * in an [optimal sorting network](https://en.wikipedia.org/wiki/Sorting_network)
 * of appropriate size.
 *
 * | Geometry type | number of vertices | min. index comparisons | max. index comparisons | optimal sorting network |
 * |---------------|--------------------|------------------------|------------------------|-------------------------|
 * | line          | 2                  | 0                      | 0                      | 1                       |
 * | triangle      | 3                  | 3                      | 3                      | 3                       |
 * | quadrilateral | 4                  | 5                      | 6                      | 5                       |
 * | tetrahedron   | 4                  | 6                      | 6                      | 5                       |
 * | pyramid       | 5                  | 9                      | 10                     | 9                       |
 * | prism         | 6                  | >=12                   | <=15                   | 12                      |
 * | hexahedron    | 8                  | >=18                   | <=24                   | 19                      |
 *
 */
template<std::size_t dim>
class FaceOrientations
{
  static constexpr std::size_t maxEdges = 12;
  static constexpr std::size_t maxFacets = 6;
  static constexpr std::size_t bitsPerFacet = 4;

  using Data = std::bitset<maxEdges + maxFacets*bitsPerFacet>;

  static constexpr Data singleFacetMask = (1<<bitsPerFacet) - 1;

  // Offset of the 2d facet within the bitfield
  static constexpr auto facetBitOffset(int subEntity)
  {
    return maxEdges + bitsPerFacet*subEntity;
  }

  template<unsigned int blockSize, class T>
  static constexpr T getBitBlock(const T& bits, unsigned int i)
  {
    constexpr T blockMask = (1 << (blockSize)) - 1;
    return (bits >> (blockSize*i)) & blockMask;
  }

  /**
   * \brief Compute orientation of line-face with given index from range global vertex indices
   *
   * \param vertexIndices Range of global indices of the face
   * \param re Reference element of the entity
   * \param subEntity Local index of the line-face
   *
   * This computes the line orientation directly from the vertex indices.
   */
  constexpr void computeEdgeOrientation(const auto& vertexIndices, const auto& re, int subEntity)
  {
    data_[subEntity] = vertexIndices[0] > vertexIndices[1];
  }

  /**
   * \brief Compute orientation of triangle-face with given index from range global vertex indices
   *
   * \param vertexIndices Range of global indices of the face
   * \param re Reference element of the entity
   * \param subEntity Local index of the triangle-face
   *
   * This assumes the orientations of the edges of the triangle have already
   * been computed and uses these orientations to avoid index comparisons.
   * Since the triangle orientation is fully determined by the edge orientations,
   * no further index comparisons are required.
   */
  constexpr void computeTriangleOrientation(const auto& vertexIndices, const auto& re, int subEntity)
  {
    // All possible triangle orientations are stored in a binary encoded lookup table.
    // We map the orientations of the three adjacent edges to a flat
    // index in 0,...,7 and use it to access a 4 bit block in the
    // table encoding the triangle orientation.
    // The leading bit is always zero to indicate that the face is a triangle.
    constexpr uint32_t edgeOrientationsToTriangleOrientation = 0b0101'0001'0011'0100'0010'0011'0110'0000;
    auto edgeOrientations = (data_[re.subEntity(subEntity, 1, 0, 2)] << 0)
                          | (data_[re.subEntity(subEntity, 1, 1, 2)] << 1)
                          | (data_[re.subEntity(subEntity, 1, 2, 2)] << 2);
    auto faceOrientation = getBitBlock<4>(edgeOrientationsToTriangleOrientation, edgeOrientations);
    // Set bits associated to this facet.
    data_ |= (Data(faceOrientation) << facetBitOffset(subEntity));
  }

  /**
   * \brief Compute orientation of quadrilateral-face with given index from range global vertex indices
   *
   * \param vertexIndices Range of global indices of the face
   * \param subEntity Local index of the quadrilateral-face
   *
   * This assumes the orientations of the edges of the quadrilateral have already
   * been computed and uses these orientations to avoid index comparisons.
   * Since the quadrilateral orientation is not fully determined by the edge orientations,
   * this requires additional index comparisons.
   */
  constexpr void computeQuadrilateralOrientation(const auto& vertexIndices, const auto& re, auto subEntity)
  {
    // The following computes the local index i_min of the quadrilateral
    // vertex with smallest index. I.e. it is equivalent to
    //
    //   std::size_t i_min = 0;
    //   for(auto i: Dune::range(1, 4))
    //     if (vertexIndices[i] < vertexIndices[i_min])
    //       i_min = i;
    //
    // This code would always do 3 index comparisons although several of them have
    // already been done for the edges. Hence we rather want to reuse the latter
    // to save comparisons (for large index/id types). The following code maps
    // the 16 possible outcomes of the edge comparisons to the index i which
    // is then used to lookup i_min from a precomputed index table. Since i_min
    // is from 0,..,3 we can represent the value by two bits and encode the
    // lookup table as 16*2 bits. There are two cases (5 and 10), where i_min
    // cannot be deduced from edge comparisons because the two vertices with the
    // smallest index are not neighbors. These cases require an additional comparison
    // of the indices of those two vertices.
    //
    // In contrast to the linear search above this reduces the total number of comparisons
    // for a hexahedron from 36=12+6*4 to a number between 18=12+6 and 24=12+6*2,
    // depending on the number of required non-neighbor comparisons.
    //
    // For index/id types that are cheap to compare like e.g. int, the code does
    // not make a measurable difference. For larger types it is slightly faster.
    constexpr uint32_t edgeOrientationsToMinVertex = 0b11'11'01'01'11'00'00'00'10'00'00'01'10'00'10'00;
    auto edgeOrientations = (data_[re.subEntity(subEntity, 1, 0, 2)] << 0)
                          | (data_[re.subEntity(subEntity, 1, 1, 2)] << 1)
                          | (data_[re.subEntity(subEntity, 1, 2, 2)] << 2)
                          | (data_[re.subEntity(subEntity, 1, 3, 2)] << 3);
    std::size_t i_min = 0;
    // In case 5 and 10 we need to do an extra comparison, because
    // we cannot deduce i_min from edge orientations. In all other
    // cases, we can lookup i_min by extracting two bits from the table.
    if (edgeOrientations==5)
      i_min = (vertexIndices[1] < vertexIndices[2]) ? 1 : 2;
    else if (edgeOrientations==10)
      i_min = (vertexIndices[0] < vertexIndices[3]) ? 0 : 3;
    else
      i_min = getBitBlock<2>(edgeOrientationsToMinVertex, edgeOrientations);
    // The lowest two bits of the faceOrientation indicate the number of
    // counter clockwise rotations needed to map vertex i_min to vertex 0.
    // The third bit indicates if we need a reflection to revert the order.
    // The leading bit is always one to indicate that the face is a quadrilateral.
    unsigned long faceOrientation = 0;
    if (i_min==0)
      faceOrientation = 0b1000 | ((vertexIndices[2] < vertexIndices[1]) << 2);
    else if (i_min==1)
      faceOrientation = 0b1011 | ((vertexIndices[0] < vertexIndices[3]) << 2);
    else if (i_min==2)
      faceOrientation = 0b1001 | ((vertexIndices[3] < vertexIndices[0]) << 2);
    else if (i_min==3)
      faceOrientation = 0b1010 | ((vertexIndices[1] < vertexIndices[2]) << 2);
    // Set bits associated to this facet.
    data_ |= (Data(faceOrientation) << facetBitOffset(subEntity));
  }

public:

  /**
   * \brief Default construct FaceOrientations with index 0
   */
  FaceOrientations() = default;

  /**
   * \brief Construct FaceOrientations for given vertex ids
   *
   * \param geometryType GeometryType of the element
   * \param vertexIds Range of global vertex ids for all element vertices
   * \param codimensions Dune::index_constants for all requested codimensions
   *
   * The orientation will only be computed for sub-entities of the requested
   * codimensions. Codimensions larger than the grid dimension will be ignored.
   */
  template <class VertexIds, std::size_t... codims>
  FaceOrientations(const Dune::GeometryType& geometryType, const VertexIds& vertexIds, const Dune::index_constant<codims>&... codimensions)
    : data_(0)
  {
    if constexpr ((dim == 1) or (sizeof...(codims)==0))
      return;

    // Get reference element
    const auto& re = Dune::referenceElement<double,dim>(geometryType);

    // Access and store vertex ids only once
    using Index = std::decay_t<decltype(vertexIds[0])>;
    auto cachedVertexIds = std::array<Index, Dune::power(2, dim)>{};
    for(auto k : Dune::range(re.size(dim)))
      cachedVertexIds[k] = vertexIds[k];

    // Subrange of face vertex ids
    auto faceVertexIds = [&](auto faceIndex, auto faceCodim) {
      return Dune::transformedRangeView(re.subEntities(faceIndex, faceCodim, dim), [&](auto localVertexIndex) {
        return cachedVertexIds[localVertexIndex];
      });
    };

    if constexpr ((dim > 1) and ((codims<=(dim-1)) || ...))
    {
      for(auto subEntity : Dune::range(re.size(dim-1)))
        computeEdgeOrientation(faceVertexIds(subEntity, dim-1), re, subEntity);
    }
    if constexpr ((dim == 3) and ((codims==(dim-2)) || ...))
    {
      for(auto subEntity : Dune::range(re.size(1)))
      {
        if (re.type(subEntity, 1).isTriangle())
          computeTriangleOrientation(faceVertexIds(subEntity, 1), re, subEntity);
        if (re.type(subEntity, 1).isQuadrilateral())
          computeQuadrilateralOrientation(faceVertexIds(subEntity, 1), re, subEntity);
      }
    }
  }

  /**
   * \brief Return index of the orientation for the given face
   */
  std::size_t faceOrientationIndex(std::size_t subEntity, std::size_t codim) const
  {
    if (codim == dim-1)
      return data_[subEntity];
    if (codim == 1)
      return ((data_ >> facetBitOffset(subEntity)) & singleFacetMask).to_ulong();
    return 0;
  }

private:

  Data data_;
};

} // namespace Experimental



// *****************************************************************************
// This is the reusable part of the LagrangeBasis. It contains
//
//   LagrangePreBasis
//   LagrangeNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename R=double>
class LagrangeNode;

template<typename GV, int k, typename R=double>
class LagrangePreBasis;



/**
 * \brief A pre-basis for a PQ-lagrange bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   Range type used for shape function values
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \warning For pyramid elements in 3d, the shape functions are different for
 *    run-time and compile-time order. While the former are defined using the
 *    Duffy-transformation the latter are continuous and piecewise polynomial
 *    with discontinuous gradients.
 */
template<typename GV, int k, typename R>
class LagrangePreBasis :
  public LeafPreBasisMixin< LagrangePreBasis<GV,k,R> >
{
  static const int dim = GV::dimension;
  static const bool useDynamicOrder = (k<0);

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = LagrangeNode<GV, k, R>;

  //! Constructor for a given grid view object with compile-time order
  LagrangePreBasis(const GridView& gv)
  : LagrangePreBasis(gv, std::numeric_limits<unsigned int>::max())
  {}

  //! Constructor for a given grid view object and run-time order
  LagrangePreBasis(const GridView& gv, unsigned int order) :
    gridView_(gv), order_(order)
  {
    if (!useDynamicOrder && order!=std::numeric_limits<unsigned int>::max())
      DUNE_THROW(RangeError, "Template argument k has to be -1 when supplying a run-time order!");

    for (int i=0; i<=dim; i++)
    {
      dofsPerCube_[i] = computeDofsPerCube(i);
      dofsPerSimplex_[i] = computeDofsPerSimplex(i);
    }
    dofsPerPrism_ = computeDofsPerPrism();
    dofsPerPyramid_ = computeDofsPerPyramid();
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    vertexOffset_        = 0;
    edgeOffset_            = vertexOffset_          + dofsPerCube(0) * ((size_type)gridView_.size(dim));

    if (dim>=2)
    {
      triangleOffset_      = edgeOffset_            + dofsPerCube(1) * ((size_type) gridView_.size(dim-1));

      quadrilateralOffset_ = triangleOffset_        + dofsPerSimplex(2) * ((size_type)gridView_.size(Dune::GeometryTypes::triangle));
    }

    if (dim==3) {
      tetrahedronOffset_   = quadrilateralOffset_ + dofsPerCube(2) * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral));

      prismOffset_         = tetrahedronOffset_   +   dofsPerSimplex(3) * ((size_type)gridView_.size(Dune::GeometryTypes::tetrahedron));

      hexahedronOffset_    = prismOffset_         +   dofsPerPrism() * ((size_type)gridView_.size(Dune::GeometryTypes::prism));

      pyramidOffset_       = hexahedronOffset_    +   dofsPerCube(3) * ((size_type)gridView_.size(Dune::GeometryTypes::hexahedron));
    }
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{order()};
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    switch (dim)
    {
      case 1:
        return dofsPerCube(0) * ((size_type)gridView_.size(dim))
          + dofsPerCube(1) * ((size_type)gridView_.size(dim-1));
      case 2:
      {
        return dofsPerCube(0) * ((size_type)gridView_.size(dim))
          + dofsPerCube(1) * ((size_type)gridView_.size(dim-1))
          + dofsPerSimplex(2) * ((size_type)gridView_.size(Dune::GeometryTypes::triangle))
          + dofsPerCube(2) * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral));
      }
      case 3:
      {
        return dofsPerCube(0) * ((size_type)gridView_.size(dim))
          + dofsPerCube(1) * ((size_type)gridView_.size(dim-1))
          + dofsPerSimplex(2) * ((size_type)gridView_.size(Dune::GeometryTypes::triangle))
          + dofsPerCube(2) * ((size_type)gridView_.size(Dune::GeometryTypes::quadrilateral))
          + dofsPerSimplex(3) * ((size_type)gridView_.size(Dune::GeometryTypes::tetrahedron))
          + dofsPerPyramid() * ((size_type)gridView_.size(Dune::GeometryTypes::pyramid))
          + dofsPerPrism() * ((size_type)gridView_.size(Dune::GeometryTypes::prism))
          + dofsPerCube(3) * ((size_type)gridView_.size(Dune::GeometryTypes::hexahedron));
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    // That cast to unsigned int is necessary because GV::dimension is an enum,
    // which is not recognized by the power method as an integer type...
    return power(order()+1, (unsigned int)GV::dimension);
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    for (size_type i = 0, end = node.finiteElement().size() ; i < end ; ++it, ++i)
    {
      Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
      const auto& gridIndexSet = gridView().indexSet();
      const auto& element = node.element();

      // The dimension of the entity that the current dof is related to
      auto dofDim = dim - localKey.codim();

      // Test for a vertex dof
      // The test for k==1 is redundant, but having it here allows the compiler to conclude
      // at compile-time that the dofDim==0 case is the only one that will ever happen.
      // This leads to measurable speed-up: see
      //   https://gitlab.dune-project.org/staging/dune-functions/issues/30
      if (k==1 || dofDim==0) {
        *it = {{ (size_type)(gridIndexSet.subIndex(element,localKey.subEntity(),dim)) }};
        continue;
      }

      if (dofDim==1)
        {  // edge dof
          if (dim==1)  // element dof -- any local numbering is fine
            {
              *it = {{ edgeOffset_
                       + dofsPerCube(1) * ((size_type)gridIndexSet.subIndex(element,0,0))
                       + localKey.index() }};
              continue;
            }
          else
            {
              const auto refElement
                = Dune::referenceElement<double,dim>(element.type());

              // We have to reverse the numbering if the local element edge is
              // not aligned with the global edge.
              auto v0 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
              auto v1 = (size_type)gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
              bool flip = (v0 > v1);
              *it = {{ (flip)
                       ? edgeOffset_
                       + dofsPerCube(1)*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                       + (dofsPerCube(1)-1)-localKey.index()
                       : edgeOffset_
                       + dofsPerCube(1)*((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                       + localKey.index() }};
              continue;
            }
        }

      if (dofDim==2)
        {
          if (dim==2)   // element dof -- any local numbering is fine
            {
              if (element.type().isTriangle())
                {
                  *it = {{ triangleOffset_ + dofsPerSimplex(2)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else if (element.type().isQuadrilateral())
                {
                  *it = {{ quadrilateralOffset_ + dofsPerCube(2)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else
                DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
            } else
            {
              const auto refElement
                = Dune::referenceElement<double,dim>(element.type());

              if (order()>3)
                DUNE_THROW(Dune::NotImplemented, "LagrangeBasis for 3D grids is only implemented if k<=3");

              if (order()==3 and !refElement.type(localKey.subEntity(), localKey.codim()).isTriangle())
                DUNE_THROW(Dune::NotImplemented, "LagrangeBasis for 3D grids with k==3 is only implemented if the grid is a simplex grid");

              *it = {{ triangleOffset_ + ((size_type)gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())) }};
              continue;
            }
        }

      if (dofDim==3)
        {
          if (dim==3)   // element dof -- any local numbering is fine
            {
              if (element.type().isTetrahedron())
                {
                  *it = {{ tetrahedronOffset_ + dofsPerSimplex(3)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else if (element.type().isHexahedron())
                {
                  *it = {{ hexahedronOffset_ + dofsPerCube(3)*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else if (element.type().isPrism())
                {
                  *it = {{ prismOffset_ + dofsPerPrism()*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else if (element.type().isPyramid())
                {
                  *it = {{ pyramidOffset_ + dofsPerPyramid()*((size_type)gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                  continue;
                }
              else
                DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedra, hexahedra, prisms, or pyramids");
            } else
            DUNE_THROW(Dune::NotImplemented, "Grids of dimension larger than 3 are no supported");
        }
      DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the LagrangeBasis");
    }
    return it;
  }

  //! Polynomial order used in the local Lagrange finite-elements
  unsigned int order() const
  {
    return (useDynamicOrder) ? order_ : k;
  }

protected:
  GridView gridView_;

  // Run-time order, only valid if k<0
  unsigned int order_;

  //! Number of degrees of freedom assigned to a simplex (without the ones assigned to its faces!)
  size_type dofsPerSimplex(std::size_t simplexDim) const
  {
    return useDynamicOrder ? dofsPerSimplex_[simplexDim] : computeDofsPerSimplex(simplexDim);
  }

  //! Number of degrees of freedom assigned to a cube (without the ones assigned to its faces!)
  size_type dofsPerCube(std::size_t cubeDim) const
  {
    return useDynamicOrder ? dofsPerCube_[cubeDim] : computeDofsPerCube(cubeDim);
  }

  size_type dofsPerPrism() const
  {
    return useDynamicOrder ? dofsPerPrism_ : computeDofsPerPrism();
  }

  size_type dofsPerPyramid() const
  {
    return useDynamicOrder ? dofsPerPyramid_ : computeDofsPerPyramid();
  }

  //! Number of degrees of freedom assigned to a simplex (without the ones assigned to its faces!)
  size_type computeDofsPerSimplex(std::size_t simplexDim) const
  {
    return order() == 0 ? (dim == simplexDim ? 1 : 0) : Dune::binomial(std::size_t(order()-1),simplexDim);
  }

  //! Number of degrees of freedom assigned to a cube (without the ones assigned to its faces!)
  size_type computeDofsPerCube(std::size_t cubeDim) const
  {
    return order() == 0 ? (dim == cubeDim ? 1 : 0) : Dune::power(order()-1, cubeDim);
  }

  size_type computeDofsPerPrism() const
  {
    return order() == 0 ? (dim == 3 ? 1 : 0) : (order()-1)*(order()-1)*(order()-2)/2;
  }

  size_type computeDofsPerPyramid() const
  {
    return order() == 0 ? (dim == 3 ? 1 : 0) : (order()-2)*(order()-1)*(2*order()-3)/6;
  }

  // When the order is given at run-time, the following numbers are pre-computed:
  std::array<size_type,dim+1> dofsPerSimplex_;
  std::array<size_type,dim+1> dofsPerCube_;
  size_type dofsPerPrism_;
  size_type dofsPerPyramid_;

  size_type vertexOffset_;
  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type quadrilateralOffset_;
  size_type tetrahedronOffset_;
  size_type pyramidOffset_;
  size_type prismOffset_;
  size_type hexahedronOffset_;

};



template<typename GV, int k, typename R>
class LagrangeNode :
  public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

  // A simple cache handing storing exactly one LFE.  This can be
  // used for grids with only a single element type. In contrast to
  // StaticLagrangeLocalFiniteElementCache this also supports
  // Lagrange*LocalFiniteElement with run-time order.
  template <class LFE>
  class SingleLocalFiniteElementCache
  {
    LFE lfe_;
  public:
    using FiniteElementType = LFE;

    template<class... Args>
    SingleLocalFiniteElementCache(Args&&... args)
      : lfe_(std::forward<Args>(args)...)
    {}

    //! Obtain the cached local finite-element.
    const FiniteElementType& get ([[maybe_unused]] Dune::GeometryType type) const
    {
      return lfe_;
    }
  };

  // Utility function to construct the FiniteElementCache type.
  // Since the function is just a helper to generate a type,
  // it hands out a MetaType<T> instead of a raw T.
  static constexpr auto makeCacheType()
  {
    using D = typename GV::ctype;
    if constexpr (Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>::v)
    {
      constexpr auto type = Dune::GeometryType(Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId, GV::dimension);
      if constexpr (type.isSimplex())
        return Dune::MetaType<SingleLocalFiniteElementCache<Dune::LagrangeSimplexLocalFiniteElement<D,R,dim,k>>>{};
      else if constexpr (type.isCube())
        return Dune::MetaType<SingleLocalFiniteElementCache<Dune::LagrangeCubeLocalFiniteElement<D,R,dim,k>>>{};
      else if constexpr (type.isPrism())
        return Dune::MetaType<SingleLocalFiniteElementCache<Dune::LagrangePrismLocalFiniteElement<D,R,k>>>{};
      else if constexpr (type.isPyramid())
        return Dune::MetaType<SingleLocalFiniteElementCache<Dune::LagrangePyramidLocalFiniteElement<D,R,k>>>{};
    }
    else
      return Dune::MetaType<Dune::LagrangeLocalFiniteElementCache<D,R,dim,k>>{};
  }

  using FiniteElementCache = typename decltype(makeCacheType())::type;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  //! Constructor without order (uses the compile-time value)
  LagrangeNode() :
    LagrangeNode(k)
  {}

  //! Constructor with a run-time order
  LagrangeNode(unsigned int order) :
    cache_(order),
    finiteElement_(nullptr),
    element_(nullptr)
  {}

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
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
    this->setSize(finiteElement_->size());
  }

protected:

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   The range type of the local basis
 */
template<std::size_t k, typename R=double>
auto lagrange()
{
  return [](const auto& gridView) {
    return LagrangePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
}

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis with a run-time order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 */
template<typename R=double>
auto lagrange(int order)
{
  return [=](const auto& gridView) {
    return LagrangePreBasis<std::decay_t<decltype(gridView)>, -1, R>(gridView, order);
  };
}

} // end namespace BasisFactory



/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of LagrangePreBasis.
 *
 * \warning The implementation of the basis with run-time order order uses the
 *   LagrangeFiniteElement implementation of dune-localfunctions, which is known
 *   to violate strict-aliasing rules
 *   (see https://gitlab.dune-project.org/core/dune-localfunctions/issues/14)
 *   Keep this in mind if ever you experience difficult-to-explain crashes
 *   or wrong results.
 *
 * \warning For pyramid elements in 3d, the shape functions are different for
 *    run-time and compile-time order. While the former are defined using the
 *    Duffy-transformation the latter are continuous and piecewise polynomial
 *    with discontinuous gradients.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis; -1 means 'order determined at run-time'
 * \tparam R The range type of the local basis
 */
template<typename GV, int k=-1, typename R=double>
using LagrangeBasis = DefaultGlobalBasis<LagrangePreBasis<GV, k, R> >;





} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
