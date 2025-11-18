// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARGYRISBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARGYRISBASIS_HH

#include <algorithm>
#include <array>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/functions/analyticfunctions/monomialset.hh>
#include <dune/functions/common/mapperutilities.hh>
#include <dune/functions/functionspacebases/cubichermitebasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/functionaldescriptor.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/transformedfiniteelementmixin.hh>

/**
 * \file Argyrisbasis.hh
 * \brief This file provides an implementation of the cubic Argyris finite element.
 *
 * For reference, see[Ciarlet, The Finite Element Method for Elliptic Problems, 2002].
 * It contains in the following order:
 *     - A GlobalBasis typedef ArgyrisBasis
 *     - A template ArgyrisLocalFiniteElement providing an implementation
 *       of the LocalFiniteElement interface, along with its subparts (Impl namespace)
 *     - A template ArgyrisNode
 *     - A template ArgyrisPreBasis
 *     - A Factory argyris() in the BasisFactory namespace
 */
namespace Dune::Functions
{

  template<class GV, class R>
  class ArgyrisPreBasis;

  /** \brief Nodal basis of a scalar cubic Argyris finite element space
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \note This only works for simplex grids. The Argyris basis is only implemented for degree 5.
   * \note The Argyris Finite element has the following properties:
   *   - Its global space is in \f$ C^1 \f$.
   *   - Its interpolation evaluates derivatives up to order 2, i.e. you cannot interpolate into a lambda function.
   *   - Strongly enforcing boundary conditions is not as simple as with langrange bases
   *   - It global space is not nested, i.e. the space on a refined grid is not a subspace of the
   *     space on the coarser grid.
   * All arguments passed to the constructor will be forwarded to the constructor
   * of ArgyrisPreBasis.
   *
   * \tparam GV The GridView that the space is defined on
   * \tparam R The range type of the local basis
   */
  template<class GV, class R = double>
  using ArgyrisBasis = DefaultGlobalBasis<ArgyrisPreBasis<GV, R> >;

  namespace Impl
  {

    /** \brief Associations of the Argyris degrees of freedom to subentities of the triangle
     */
    class ArgyrisLocalCoefficients
    {
    public:
      using size_type = std::size_t;
      static constexpr int dim = 2;

      ArgyrisLocalCoefficients()
      {
        for (unsigned int i = 0; i < 3; i++)             // subentities: three vertices
          for (unsigned int k = 0; k < 6; k++)           // 6 basis functions per vertex
            localKeys_[6 * i + k]
              = LocalKey(i, dim, k);   //(subentity, codim, number of dof for this subentity)
        for (unsigned int i = 0; i < 3; ++i)             // subentities: three edges
          localKeys_[18 + i] = LocalKey(i, dim - 1, 0);  // one node per edge
      }

      /** \brief number of coefficients
       */
      static constexpr size_type size()
      {
        return 21;
      }

      /** \brief get i'th index
       */
      LocalKey const &localKey(size_type i) const
      {
        return localKeys_[i];
      }

    private:
      std::array<Dune::LocalKey, 21> localKeys_;
    };


    /** \brief Implementation of Argyris Polynomials
     * \tparam D Type to represent the field in the domain
     * \tparam R Type to represent the field in the range
     * \tparam dim Dimension of the domain simplex
     */
    template<class D, class R>
    class ArgyrisReferenceLocalBasis
    {
      static constexpr int dim = 2;
    public:
      using Traits = H2LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>,
          FieldMatrix<R, 1, dim>, FieldMatrix<R, dim, dim> >;
    private:

      /**
       * \brief Get the Argyris Coefficients Matrix
       * \return FieldMatrix<F, size, size>
       *  where size is the dimension of the quintic polynomial space (21)
       *
       * This returns the basis transformation matrix from a monomial
       * basis to the Argyris basis on the reference domain.
       * I.e. the i-th row of the returned matrix contains the
       * coefficients of the i-th Argyris basis function with respect to a
       * monomial basis. The monomials are enumerated as in the
       * MonomialSet of the corresponding dimension and order.
       * The values are taken (and reordered) from https://defelement.com/elements/argyris.html.

       */
      static constexpr auto getArgyrisCoefficients()
      {
        // Define std::sqrt(2.) manually in double precision,
        // because std::sqrt is not constexpr before C++26.
        D sqrt2 = -8. * 1.414213562373095;
        return Dune::FieldMatrix<D, 21,21>{
          // vertex functionals
          // l_0
          {1,                          /*0th order*/
           0, 0,                       /*1th order*/
           0, 0, 0,                    /*2th order*/
           -10, 0, 0, -10,             /*3th order*/
           15, 0, -30, 0, 15,          /*4th order*/
           -6, 0, 30, 30, 0, -6},      /*5th order*/
          // l_1
          {0,
           1, 0,
           0, 0, 0,
           -6, 0, -11, 0,
           8, 0, 10, 18, 0,
           -3, 0, 1, -10, -8, 0},
          // l_2, l_1 mirrored
          {
            0,
            0, 1,
            0, 0, 0,
            0, -11, 0, -6,
            0, 18, 10, 0, 8,
            0, -8, -10, 1, 0, -3},
          // l_3
          {
            0,
            0, 0,
            0.5, 0, 0,
            -1.5, 0, 0, 0,
            1.5, 0, -1.5, 0, 0,
            -0.5, 0, 1.5, 1, 0, 0},
          // l_4
          {
            0,
            0, 0,
            0, 1, 0,
            0, -4, -4, 0,
            0, 5, 10, 5, 0,
            0, -2, -6, -6, -2, 0},
          // l_5, l_3 mirrored
          {
            0,
            0, 0,
            0, 0, 0.5,
            0, 0, 0, -1.5,
            0, 0, -1.5, 0, 1.5,
            0, 0, 1, 1.5, 0, -0.5},
          // l_6
          {0,
           0, 0,
           0, 0, 0,
           10, 0, 0, 0,
           -15, 0, 15, 0, 0,
           6, 0, -15, -15, 0, 0},
          // l_7
          {0,
           0, 0,
           0, 0, 0,
           -4, 0, 0, 0,
           7, 0, -3.5, 0, 0,
           -3, 0, 3.5, 3.5, 0, 0},
          // l_8
          {0,
           0, 0,
           0, 0, 0,
           0, -5, 0, 0,
           0, 14, 18.5, 0, 0,
           0, -8, -18.5, -13.5, 0, 0},
          // l_9
          {
            0,
            0, 0,
            0, 0, 0,
            0.5, 0, 0, 0,
            -1, 0, 0.25, 0, 0,
            0.5, 0, -0.25, -0.25, 0, 0},
          // l_10
          {
            0,
            0, 0,
            0, 0, 0,
            0, 1, 0, 0,
            0, -3, -3.5, 0, 0,
            0, 2, 3.5, 2.5, 0, 0},
          // l_11
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, 0, 0,
            0, 0, 1.25, 0, 0,
            0, 0, -0.75, -1.25, 0, 0},
          // l_12 mirrors l_6
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, 0, 10,
            0, 0, 15, 0, -15,
            0, 0, -15, -15, 0, 6},
          // l_13 mirrors l_8
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, -5, 0,
            0, 0, 18.5, 14, 0,
            0, 0, -13.5, -18.5, -8, 0},
          // l_14 mirrors l_7
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, 0, -4,
            0, 0, -3.5, 0, 7,
            0, 0, 3.5, 3.5, 0, -3},
          // l_15 mirrors l_11
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, 0, 0,
            0, 0, 1.25, 0, 0,
            0, 0, -1.25, -0.75, 0, 0},
          // l_16 mirrors l_10
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, 1, 0,
            0, 0, -3.5, -3, 0,
            0, 0, 2.5, 3.5, 2, 0},
          // l_17 mirrors l_9
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, 0, 0.5,
            0, 0, 0.25, 0, -1,
            0, 0, -0.25, -0.25, 0, 0.5},
          // edge functionals
          // l_18
          {
            0,
            0, 0,
            0, 0, 0,
            0, 16, 0, 0,
            0, -32, -32, 0, 0,
            0, 16, 32, 16, 0, 0},
          // l_19
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, -16, 0,
            0, 0, 32, 32, 0,
            0, 0, -16, -32, -16, 0},
          // l_20
          {
            0,
            0, 0,
            0, 0, 0,
            0, 0, 0, 0,
            0, 0, -1. * sqrt2, 0, 0,
            0, 0, sqrt2, sqrt2, 0, 0}};
      }


      static constexpr auto referenceBasisCoefficients = getArgyrisCoefficients();
      static constexpr MonomialSet<typename Traits::RangeFieldType, dim, 5> monomials = {};

    public:

      /** The number of basis functions in the basis
       */
      static constexpr unsigned int size()
      {
        return ArgyrisLocalCoefficients::size();
      }

      /** The polynomial order of the basis
       */
      static constexpr unsigned int order()
      {
        return 5;
      }

      /** \brief Evaluate function values of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Values of all shape functions at that point
       */
      void evaluateFunction(const typename Traits::DomainType &in,
                            std::vector<typename Traits::RangeType> &out) const
      {
        out.resize(size());
        auto monomialValues = monomials(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
      }

      /** \brief Evaluate Jacobians of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Jacobians of all shape functions at that point
       */
      void evaluateJacobian(const typename Traits::DomainType &in,
                            std::vector<typename Traits::JacobianType> &out) const
      {
        out.resize(size());
        auto monomialValues = derivative(monomials)(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
      }

      /** \brief Evaluate Hessians of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Hessians of all shape functions at that point
       */
      void evaluateHessian(const typename Traits::DomainType &in,
                           std::vector<typename Traits::HessianType> &out) const
      {
        out.resize(size());
        auto monomialValues = derivative(derivative(monomials))(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
      }

      /** \brief Evaluate partial derivatives of all shape functions at a given point
       *
       * \param[in] order The partial derivative to be computed, as a multi-index
       * \param[in] in  The evaluation point
       * \param[out] out Jacobians of all shape functions at that point
       */
      void partial(std::array<unsigned int, dim> order, const typename Traits::DomainType &in,
                   std::vector<typename Traits::RangeType> &out) const
      {
        out.resize(size());
        auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
        if (totalOrder == 0)
          evaluateFunction(in, out);
        else if (totalOrder == 1)
        {
          evaluateJacobian(in,jacobiansBuffer_);
          std::size_t which = std::max_element(order.begin(), order.end()) - order.begin();
          for (auto i : Dune::range(size()))
            out[i] = jacobiansBuffer_[i][0][which];
        }
        else if (totalOrder == 2)
        {
          evaluateHessian(in, hessianBuffer_);
          std::size_t first, second;
          first = std::max_element(order.begin(), order.end()) - order.begin();
          if (order[first] == 2)
            second = first;
          else
          {
            order[first] = 0;
            second = std::max_element(order.begin(), order.end()) - order.begin();
          }
          for (auto i : Dune::range(size()))
            out[i] = hessianBuffer_[i][first][second];
        }
        else
          DUNE_THROW(RangeError, "partial() not implemented for given order");
      }

    private:
      mutable std::vector<typename Traits::JacobianType> jacobiansBuffer_;
      mutable std::vector<typename Traits::HessianType> hessianBuffer_;
    };


    /** \brief Class that evaluates the push forwards of the global nodes of a
     * LocalFunction. It stretches the LocalInterpolation interface, because we
     * evaluate the derivatives of f.
     *
     */
    template<class D>
    class ArgyrisLocalInterpolation
    {
      using size_type = std::size_t;
      static constexpr int dim = 2;

      static constexpr unsigned int size()
      {
        return ArgyrisLocalCoefficients::size();
      }

      using FunctionalDescriptor = Dune::Functions::Impl::FunctionalDescriptor<dim>;

    public:

      ArgyrisLocalInterpolation()
      {
        // first vertex
        // function evaluation
        descriptors_[0] = FunctionalDescriptor();
        // gradient evaluations
        descriptors_[1] = FunctionalDescriptor({1,0});
        descriptors_[2] = FunctionalDescriptor({0,1});
        // hessian evaluations
        descriptors_[3] = FunctionalDescriptor({2,0});
        descriptors_[4] = FunctionalDescriptor({1,1});
        descriptors_[5] = FunctionalDescriptor({0,2});
        // second vertex
        descriptors_[6] = FunctionalDescriptor();
        descriptors_[7] = FunctionalDescriptor({1,0});
        descriptors_[8] = FunctionalDescriptor({0,1});
        descriptors_[9]  = FunctionalDescriptor({2,0});
        descriptors_[10] = FunctionalDescriptor({1,1});
        descriptors_[11] = FunctionalDescriptor({0,2});
        // third vertex
        descriptors_[12] = FunctionalDescriptor();
        descriptors_[13] = FunctionalDescriptor({1,0});
        descriptors_[14] = FunctionalDescriptor({0,1});
        descriptors_[15] = FunctionalDescriptor({2,0});
        descriptors_[16] = FunctionalDescriptor({1,1});
        descriptors_[17] = FunctionalDescriptor({0,2});
        // normal derivatives at edge midpoints
        descriptors_[18] = FunctionalDescriptor(1);
        descriptors_[19] = FunctionalDescriptor(1);
        descriptors_[20] = FunctionalDescriptor(1);
      }

      /** \brief bind the Interpolation to an element and a local state.
       */
      template<class Element>
      void bind( Element const& element, std::array<D, 3>const& averageVertexMeshSize, std::bitset<3>const& edgeOrientation)
      {
        averageVertexMeshSize_ = &averageVertexMeshSize;

        auto geometry = element.geometry();

        // get global Normals and midpoints
        auto refElement = Dune::referenceElement<D, 2>(geometry.type());
        for (std::size_t i = 0; i < 3; ++i)
        {
          localVertices_[i] = refElement.position(i, 2);

          localMidpoints_[i] = refElement.position(i, 1);
          std::size_t lower = (i == 2) ? 1 : 0;
          std::size_t upper = (i == 0) ? 1 : 2;

          auto edge = geometry.global(refElement.position(upper, 2))
                      - geometry.global(refElement.position(lower, 2));
          // normalize and orient
          edge /= edge.two_norm() * (edgeOrientation[i] ? -1. : 1.);
          // Rotation by pi/2. Note that Kirby rotates by 3*pi/2
          globalNormals_[i] = {-edge[1], edge[0]};
        }
      }

      /** \brief Evaluate a given function and its derivatives at the nodes
       *
       * \tparam F Type of function to evaluate
       * \tparam C Type used for the values of the function
       * \param[in] f Function to evaluate
       * \param[out] out Array of function values
       */
      template<class F, class C>
      void interpolate(F const& f, std::vector<C>& out) const
      {
        out.resize(size());
        auto&& df = derivative(f);
        auto&& Hf = derivative(df);

        int offset = 0;
        // iterate over vertices
        for (int i = 0; i < 3; i++, offset += 6)
        {
          auto&& x = localVertices_[i];
          auto derivativeValue = squeezeTensor(df(x));
          auto hessianValue = squeezeTensor(Hf(x));
          auto && h = (*averageVertexMeshSize_)[i];

          out[offset ] = f(x);

          out[offset + 1] = derivativeValue[0] * h;
          out[offset + 2] = derivativeValue[1] * h;

          out[offset + 3] = hessianValue[0][0] * h*h;
          out[offset + 4] = hessianValue[0][1] * h*h;
          out[offset + 5] = hessianValue[1][1] * h*h;
        }
        // iterate over edges
        for (int i = 0; i < 3; i++)
          out[18 + i] = squeezeTensor(df(localMidpoints_[i])).dot(globalNormals_[i]);

      }

      /**
       * \brief Get the object that describes which type of evaluations is performed by the dual basis
       */
      const FunctionalDescriptor& functionalDescriptor(size_type i) const
      {
        return descriptors_[i];
      }

    protected:
      std::array<Dune::FieldVector<D, 2>, 3> globalNormals_;
      std::array<Dune::FieldVector<D, 2>, 3> localMidpoints_;
      std::array<Dune::FieldVector<D, 2>, 3> localVertices_;
      std::array<D, 3> const* averageVertexMeshSize_;
      std::array<FunctionalDescriptor, size()> descriptors_;
    };

    // TODO: There is a copy of this further up
    template<class D, class R>
    struct ArgyrisLocalBasisTraits
      : public H2LocalBasisTraits<D, 2, Dune::FieldVector<D,2>, R, 1,
            Dune::FieldVector<R,1>, Dune::FieldMatrix<R,1,2>, Dune::FieldMatrix<R,2,2> >
    {};

    /** \brief Argyris finite element for triangles, as defined on the reference Element
     *
     * Note that this is a non affine-equivalent finite element.
     * It requires an additional transformation to the relate reference basis
     * with the pullbacks of global basis.
     *
     * \tparam D Type used for domain coordinates
     * \tparam R Type used for function values
     */
    template<class D, class R>
    class ArgyrisLocalFiniteElement
      : public Impl::TransformedFiniteElementMixin<ArgyrisLocalFiniteElement<D,R>, ArgyrisLocalBasisTraits<D, R> >
    {
      using Base = Impl::TransformedFiniteElementMixin< ArgyrisLocalFiniteElement<D,R>, ArgyrisLocalBasisTraits<D, R> >;
      friend class Impl::TransformedLocalBasis<ArgyrisLocalFiniteElement<D,R>, ArgyrisLocalBasisTraits<D, R> >;
      static constexpr int dim = 2;

    public:
      /** \brief Export number types, dimensions, etc.
       */
      using size_type = std::size_t;
      using Traits = LocalFiniteElementTraits<
          Impl::TransformedLocalBasis<ArgyrisLocalFiniteElement<D,R>, ArgyrisLocalBasisTraits<D, R> >,
          Impl::ArgyrisLocalCoefficients,
          Impl::ArgyrisLocalInterpolation<D> >;

      /** \brief Returns the assignment of the degrees of freedom to the element
       * subentities
       */
      const typename Traits::LocalCoefficientsType &localCoefficients() const
      {
        return coefficients_;
      }

      /** \brief Returns object that evaluates degrees of freedom
       */
      const typename Traits::LocalInterpolationType &localInterpolation() const
      {
        return interpolation_;
      }

      /** \brief The reference element that the local finite element is defined on
       */
      static constexpr GeometryType type()
      {
        return GeometryTypes::simplex(dim);
      }

      /** The size of the transformed finite element.
       */
      static constexpr size_type size()
      {
        return Impl::ArgyrisLocalCoefficients::size();
      }

      /** Binds the Finite Element to an element.
       */
      template<class VertexMapper,class ElementMapper, class Element>
      void bind(VertexMapper const& vertexMapper, std::vector<D> const& globalAverageVertexMeshSize,
                ElementMapper const& elementMapper, std::vector<std::bitset<3> >const& edgeOrientations, Element const &e)
      {
        // Cache average mesh size for each vertex
        for (auto i : range(dim+1))
          averageVertexMeshSize_[i] = globalAverageVertexMeshSize[vertexMapper.subIndex(e, i, dim)];

        // Cache orientation for each mesh
        edgeOrientation_ = edgeOrientations[elementMapper.index(e)];

        // Bind LocalInterpolation to updated local state
        interpolation_.bind(e, averageVertexMeshSize_, edgeOrientation_);

        // Compute local transformation matrices for each vertex
        const auto& geometry = e.geometry();
        const auto& refElement = ReferenceElements<typename Element::Geometry::ctype, dim>::simplex();
        for (auto i : range(dim+1))
        {
          vertexJacobians_[i] = geometry.jacobian(refElement.position(i, dim));
        }

        // get geometrical information
        std::array<FieldVector<R, 2>, 3> referenceTangents;

        // get local and global Tangents
        for (std::size_t i = 0; i < 3; ++i)
        {
          std::size_t lower = (i == 2) ? 1 : 0;
          std::size_t upper = (i == 0) ? 1 : 2;
          auto edge = refElement.position(upper, 2) - refElement.position(lower, 2);

          // store normalized reference Tangent vectors
          referenceTangents[i] = edge / edge.two_norm();

          auto globalEdge = geometry.global(refElement.position(upper, 2))
                            - geometry.global(refElement.position(lower, 2));

          // store length of global tangents and normalized global tangent vectors
          l[i] = globalEdge.two_norm();
          globalTangents[i] = globalEdge / l[i];

          tau[i] = FieldVector<R, 3>{
            globalTangents[i][0] * globalTangents[i][0],
            2. * globalTangents[i][0] * globalTangents[i][1],
            globalTangents[i][1] * globalTangents[i][1]
          };

          // The following variables are named as in Kirby's paper
          // two g matrices
          auto&& referenceG = FieldMatrix<R, 2, 2>{
            {-referenceTangents[i][1], referenceTangents[i][0]},
            {referenceTangents[i][0], referenceTangents[i][1]}
          };

          auto&& globalG = FieldMatrix<R, 2, 2>{
            {-globalTangents[i][1], globalTangents[i][0]},
            {globalTangents[i][0], globalTangents[i][1]}
          };

          //
          b[i] = globalG * geometry.jacobian(refElement.position(i, 2)) * referenceG;

          theta[i][0][0] = vertexJacobians_[i][0][0] * vertexJacobians_[i][0][0];
          theta[i][0][1] = vertexJacobians_[i][0][0] * vertexJacobians_[i][0][1];
          theta[i][0][2] = vertexJacobians_[i][0][1] * vertexJacobians_[i][0][1];

          theta[i][1][0] = 2. * vertexJacobians_[i][0][0] * vertexJacobians_[i][1][0];
          theta[i][1][1] = vertexJacobians_[i][0][0] * vertexJacobians_[i][1][1] + vertexJacobians_[i][0][1] * vertexJacobians_[i][1][0];
          theta[i][1][2] = 2. * vertexJacobians_[i][0][1] * vertexJacobians_[i][1][1];

          theta[i][2][0] = vertexJacobians_[i][1][0] * vertexJacobians_[i][1][0];
          theta[i][2][1] = vertexJacobians_[i][1][0] * vertexJacobians_[i][1][1];
          theta[i][2][2] = vertexJacobians_[i][1][1] * vertexJacobians_[i][1][1];
        }
      }

    protected:

      /** \brief Returns the local basis, i.e., the set of shape functions
       */
      Impl::ArgyrisReferenceLocalBasis<D, R> const& referenceLocalBasis() const
      {
        return basis_;
      }

      /** Applies the transformation. Note that we do not distinguish for
       * Scalar/Vector/Matrix Type, but only assume the Values to be Elements of a Vectorspace.
       * We assume containers with random access iterators.
       */
      template<class InputValues, class OutputValues>
      void transform(InputValues const &inValues, OutputValues &outValues) const
      {
        using std::pow;
        assert(inValues.size() == size());
        assert(outValues.size() == inValues.size());
        // compatibility with sympy code below
        auto &[b_0, b_1, b_2] = b;
        auto &[J_0, J_1, J_2] = vertexJacobians_;
        auto &[theta_0, theta_1, theta_2] = theta;
        auto & h = averageVertexMeshSize_;
        auto & o = edgeOrientation_;
        std::array<Dune::FieldVector<R, 2>, 3> const &t = globalTangents;
        // This code is generated with sympy.
        // It a matrix free implementation of the matrix V in Kirbys paper.
        outValues[0] = -15.0/8.0*b_0[1][0]*inValues[18]/l[0] - 15.0/8.0*b_1[1][0]*inValues[19]/l[1] + inValues[0];
        outValues[1] = (J_0[0][0]*inValues[1] + J_0[0][1]*inValues[2] - 0.4375*b_0[1][0]*inValues[18]*t[0][0] - 0.4375*b_1[1][0]*inValues[19]*t[1][0])/h[0];
        outValues[2] = (J_0[1][0]*inValues[1] + J_0[1][1]*inValues[2] - 0.4375*b_0[1][0]*inValues[18]*t[0][1] - 0.4375*b_1[1][0]*inValues[19]*t[1][1])/h[0];
        outValues[3] = (-1.0/32.0*b_0[1][0]*inValues[18]*l[0]*tau[0][0] - 1.0/32.0*b_1[1][0]*inValues[19]*l[1]*tau[1][0] + inValues[3]*theta_0[0][0] + inValues[4]*theta_0[0][1] + inValues[5]*theta_0[0][2])/h[0] /h[0];
        outValues[4] = (-1.0/32.0*b_0[1][0]*inValues[18]*l[0]*tau[0][1] - 1.0/32.0*b_1[1][0]*inValues[19]*l[1]*tau[1][1] + inValues[3]*theta_0[1][0] + inValues[4]*theta_0[1][1] + inValues[5]*theta_0[1][2])/h[0] /h[0];
        outValues[5] = (-1.0/32.0*b_0[1][0]*inValues[18]*l[0]*tau[0][2] - 1.0/32.0*b_1[1][0]*inValues[19]*l[1]*tau[1][2] + inValues[3]*theta_0[2][0] + inValues[4]*theta_0[2][1] + inValues[5]*theta_0[2][2])/h[0] /h[0];
        outValues[6] = (15.0/8.0)*b_0[1][0]*inValues[18]/l[0] - 15.0/8.0*b_2[1][0]*inValues[20]/l[2] + inValues[6];
        outValues[7] = (J_1[0][0]*inValues[7] + J_1[0][1]*inValues[8] - 0.4375*b_0[1][0]*inValues[18]*t[0][0] - 0.4375*b_2[1][0]*inValues[20]*t[2][0])/h[1];
        outValues[8] = (J_1[1][0]*inValues[7] + J_1[1][1]*inValues[8] - 0.4375*b_0[1][0]*inValues[18]*t[0][1] - 0.4375*b_2[1][0]*inValues[20]*t[2][1])/h[1];
        outValues[9] = ((1.0/32.0)*b_0[1][0]*inValues[18]*l[0]*tau[0][0] - 1.0/32.0*b_2[1][0]*inValues[20]*l[2]*tau[2][0] + inValues[9]*theta_1[0][0] + inValues[10]*theta_1[0][1] + inValues[11]*theta_1[0][2])/h[1] /h[1];
        outValues[10] = ((1.0/32.0)*b_0[1][0]*inValues[18]*l[0]*tau[0][1] - 1.0/32.0*b_2[1][0]*inValues[20]*l[2]*tau[2][1] + inValues[9]*theta_1[1][0] + inValues[10]*theta_1[1][1] + inValues[11]*theta_1[1][2])/h[1] /h[1];
        outValues[11] = ((1.0/32.0)*b_0[1][0]*inValues[18]*l[0]*tau[0][2] - 1.0/32.0*b_2[1][0]*inValues[20]*l[2]*tau[2][2] + inValues[9]*theta_1[2][0] + inValues[10]*theta_1[2][1] + inValues[11]*theta_1[2][2])/h[1] /h[1];
        outValues[12] = (15.0/8.0)*b_1[1][0]*inValues[19]/l[1] + (15.0/8.0)*b_2[1][0]*inValues[20]/l[2] + inValues[12];
        outValues[13] = (J_2[0][0]*inValues[13] + J_2[0][1]*inValues[14] - 0.4375*b_1[1][0]*inValues[19]*t[1][0] - 0.4375*b_2[1][0]*inValues[20]*t[2][0])/h[2];
        outValues[14] = (J_2[1][0]*inValues[13] + J_2[1][1]*inValues[14] - 0.4375*b_1[1][0]*inValues[19]*t[1][1] - 0.4375*b_2[1][0]*inValues[20]*t[2][1])/h[2];
        outValues[15] = ((1.0/32.0)*b_1[1][0]*inValues[19]*l[1]*tau[1][0] + (1.0/32.0)*b_2[1][0]*inValues[20]*l[2]*tau[2][0] + inValues[15]*theta_2[0][0] + inValues[16]*theta_2[0][1] + inValues[17]*theta_2[0][2])/h[2] /h[2];
        outValues[16] = ((1.0/32.0)*b_1[1][0]*inValues[19]*l[1]*tau[1][1] + (1.0/32.0)*b_2[1][0]*inValues[20]*l[2]*tau[2][1] + inValues[15]*theta_2[1][0] + inValues[16]*theta_2[1][1] + inValues[17]*theta_2[1][2])/h[2] /h[2];
        outValues[17] = ((1.0/32.0)*b_1[1][0]*inValues[19]*l[1]*tau[1][2] + (1.0/32.0)*b_2[1][0]*inValues[20]*l[2]*tau[2][2] + inValues[15]*theta_2[2][0] + inValues[16]*theta_2[2][1] + inValues[17]*theta_2[2][2])/h[2] /h[2];
        outValues[18] = b_0[0][0]*inValues[18]*(o[0] ? -1 : 1);
        outValues[19] = b_1[0][0]*inValues[19]*(o[1] ? -1 : 1);
        outValues[20] = b_2[0][0]*inValues[20]*(o[2] ? -1 : 1);
      }

    private:
      typename Impl::ArgyrisReferenceLocalBasis<D, R> basis_;
      typename Traits::LocalCoefficientsType coefficients_;
      typename Traits::LocalInterpolationType interpolation_;
      // the transformation to correct the lack of affine equivalence boils down to
      // matrix without blockstructure, because the normal derivative dofs interact nontrivially
      // with the others, for details see Kirby's paper.
      // The following geometric information is needed
      // the jacobians per vertex
      std::array<Dune::FieldMatrix<R, dim, dim>, dim+1> vertexJacobians_;
      // a three by three matrix per vertex formalizing the rotation of the hessian in voigt notation
      std::array<FieldMatrix<R, 3, 3>, 3> theta;
      // the edge lengths
      FieldVector<R, 3> l;
      // the global tangent vectors (normalized)
      std::array<FieldVector<R, 2>, 3> globalTangents;
      // Geometric quantities without trivial meaning
      std::array<FieldVector<R, 3>, 3> tau;
      std::array<FieldMatrix<R, 2, 2>, 3> b;
      // Additionally, we collect some global information
      // the local state, i.e. a collection of global information restricted to this element
      std::array<D, dim+1> averageVertexMeshSize_;
      std::bitset<3> edgeOrientation_;

    };

  } // end namespace Impl



  // *****************************************************************************
  // This is the reusable part of the basis. It contains
  //
  //   ArgyrisPreBasis
  //   ArgyrisNode
  //
  // The pre-basis allows to create the others and is the owner of possible shared
  // state. These components do _not_ depend on the global basis and local view
  // and can be used without a global basis.
  // *****************************************************************************

  template<class GV, class R>
  class ArgyrisNode
    : public LeafBasisNode
  {
    using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;
    using FiniteElement = typename Impl::ArgyrisLocalFiniteElement<typename GV::ctype, R>;

    ArgyrisNode(Mapper const& vertexMapper, Mapper const& elementMapper,
                std::vector<typename GV::ctype> const& averageVertexMeshSize,
                std::vector<std::bitset<3> > const& edgeOrientation)
      : element_(nullptr)
      , vertexMapper_(&vertexMapper)
      , elementMapper_(&elementMapper)
      , averageVertexMeshSize_(&averageVertexMeshSize)
      , edgeOrientation_(&edgeOrientation)
    {}

    //! Return current element, throw if unbound
    Element const &element() const
    {
      return *element_;
    }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the
     * dune-localfunctions module
     */
    FiniteElement const &finiteElement() const
    {
      return finiteElement_;
    }

    //! Bind to element.
    void bind(Element const &e)
    {
      element_ = &e;
      finiteElement_.bind(*vertexMapper_, *averageVertexMeshSize_,
                          *elementMapper_,* edgeOrientation_,
                          *element_ );
      this->setSize(finiteElement_.size());

    }

    //! The order of the local basis.
    unsigned int order() const { return finiteElement_.localBasis().order(); }

  protected:
    FiniteElement finiteElement_;
    Element const* element_;
    Mapper const* vertexMapper_;
    Mapper const* elementMapper_;
    std::vector<typename GV::ctype> const* averageVertexMeshSize_;
    std::vector<std::bitset<3> > const* edgeOrientation_;
  };


  /**
   * \brief A pre-basis for a Argyrisbasis
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV  The grid view that the FE basis is defined on
   * \tparam R   Range type used for shape function values
   * \note This only works for simplex grids
   */
  template<class GV, class R>
  class ArgyrisPreBasis
    : public LeafPreBasisMapperMixin<GV>
  {
    using Base = LeafPreBasisMapperMixin<GV>;
    using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
    using Element = typename GV::template Codim<0>::Entity;
    using D = typename GV::ctype;
    static const std::size_t dim = GV::dimension;

    // helper methods to assign each subentity the number of dofs. Used by the LeafPreBasisMapperMixin.
    static constexpr auto ArgyrisMapperLayout(Dune::GeometryType type, int gridDim)
    {
      assert(gridDim == 2);
      if (type.isVertex())
        return 6; // one evaluation dof and two derivative dofs and three hessian dofs per vertex
      else if (type.isLine())
        return 1;
      else if (type.isTriangle())
        return 0;
      else
        DUNE_THROW(Dune::Exception, "Invalid Geometry type for Argyris Element!");

    }

  public:
    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

    //! Template mapping root tree path to type of created tree node
    using Node = ArgyrisNode<GridView, R>;

  public:

    //! Constructor for a given grid view object
    ArgyrisPreBasis(const GV &gv)
      : Base(gv, ArgyrisMapperLayout)
      , vertexMapper_({gv, mcmgVertexLayout()})
      , elementMapper_({gv, mcmgElementLayout()})
    {
      averageVertexMeshSize_ = Impl::computeAverageSubEntityMeshSize<D>(vertexMapper_);
      edgeOrientations_ = Impl::computeEdgeOrientations(elementMapper_);
    }

    //! Update the stored grid view, to be called if the grid has changed
    void update(GridView const &gv)
    {
      Base::update(gv);
      vertexMapper_.update(this->gridView());
      elementMapper_.update(this->gridView());
      averageVertexMeshSize_ = Impl::computeAverageSubEntityMeshSize<D>(vertexMapper_);
      edgeOrientations_ = Impl::computeEdgeOrientations(elementMapper_);
    }

    //! Create tree node
    Node makeNode() const
    {
      return Node{vertexMapper_, elementMapper_, averageVertexMeshSize_, edgeOrientations_};
    }

  protected:

    SubEntityMapper vertexMapper_;
    std::vector<D> averageVertexMeshSize_;
    SubEntityMapper elementMapper_;
    std::vector<std::bitset<3> > edgeOrientations_;

  }; // class ArgyrisPreBasis

  namespace BasisFactory
  {

    /**
     * \brief construct a PreBasisFactory for the quintic Argyris Finite Element
     *
     * \tparam R RangeFieldType
     * \return the PreBasisFactory
     * \relates ArgyrisBasis
     */
    template<class R = double>
    auto argyris()
    {
      return [=](auto const &gridView) {
               return ArgyrisPreBasis<std::decay_t<decltype(gridView)>, R>(gridView);
      };
    }

  } // end namespace BasisFactory


} // end namespace Dune::Functions

#endif
