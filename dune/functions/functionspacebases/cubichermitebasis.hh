// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CUBICHERMITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CUBICHERMITEBASIS_HH

#include <algorithm>
#include <array>
#include <numeric>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/functions/common/mapperutilities.hh>
#include <dune/functions/common/squeezetensor.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/functionaldescriptor.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/transformedfiniteelementmixin.hh>



#include <dune/functions/analyticfunctions/monomialset.hh>

/**
 * \file cubichermitebasis.hh
 * \brief This file provides an implementation of the cubic Hermite finite element in 1 to 3 dimensions.
 *
 * For reference, see[Ciarlet, The Finite Element Method for Elliptic Problems, 2002].
 * It contains in the following order:
 *     - A GlobalBasis typedef CubicHermiteBasis
 *     - A template H2LocalBasisTraits, extending the dune-localfunctions
 *       LocalBasisTraits by an exported Hessiantype
 *     - A template CubicHermiteLocalFiniteElement providing an implementation
 *       of the LocalFiniteElement interface, along with its subparts (Impl namespace)
 *     - A template CubicHermiteNode
 *     - A template CubicHermitePreBasis
 *     - Two factories hermite() and reducedHermite() in the BasisFactory namespace
 */
namespace Dune::Functions
{

  template<class GV, class R, bool reduced>
  class CubicHermitePreBasis;

  /** \brief Nodal basis of a scalar cubic Hermite finite element space
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \note This only works for simplex grids. The Hermite basis is only implemented for degree 3.
   * \note The Hermite Finite element has the following properties:
   *   - Its global space is in \f$ C^1 \f$ if 1d otherwise in \f$ H^1 \f$.
   *   - The reduced version (only 2d) is part of the Discrete Kirchhoff Triangle.
   *   - Its interpolation evaluates derivatives, i.e. you cannot interpolate into a lambda function.
   *   - Strongly enforcing boundary conditions is not as simple as with langrange bases
   *   - It global space is not nested, i.e. the space on a refined grid is not a subspace of the
   *     space on the coarser grid.
   * All arguments passed to the constructor will be forwarded to the constructor
   * of CubicHermitePreBasis.
   *
   * \tparam GV The GridView that the space is defined on
   * \tparam reduced whether to use the reduced Hermite finite element (only in 2d).
   * \tparam R The range type of the local basis
   */
  template<class GV,bool reduced = false, class R = double>
  using CubicHermiteBasis = DefaultGlobalBasis<CubicHermitePreBasis<GV, R, reduced> >;

  template<class DF, int n, class D, class RF, int m, class R, class J, class H>
  struct H2LocalBasisTraits
    : public LocalBasisTraits<DF, n, D, RF, m, R, J>
  {
    /** \brief Type to represent the Hessian
     *  When \f$ \hat\phi : \mbox{IR}^n \to \mbox{IR}\f$ then HessianType
     *  is an 2D-array of m x m components where entry H[i][j] contains
     *  the derivative  \f$\partial_i \partial_j \hat\phi \f$.
     */
    using HessianType = H;
  };

  namespace Impl
  {

    // *****************************************************************************
    // * Some helper functions for building polynomial bases from monomials
    // *****************************************************************************

    /** \brief Multiply the evaluations of the monomialSet (see dune/functions/analyticfunctions/monomialset.hh) with a coefficient matrix, here FieldMatrix.
     * \tparam KCoeff The Field type of the coefficient matrix.
     * \tparam sizePolynom The number of polynomials to evaluate.
     * \tparam sizeMonom The number of monomials handed as input.
     * \tparam In input vector. We only assume container access and a correct size.
     * \tparam Out Output vector. We only assume container access and a correct size.
     * \note This function may be used to turn the return types from the dune-functions::DifferentiableFunction
     * interface into those of the dune-localfunctions::LocalBasis interface.
     */
    template<class KCoeff, int sizePolynom, int sizeMonom, class In, class Out>
    void multiplyWithCoefficentMatrix(Dune::FieldMatrix<KCoeff, sizePolynom, sizeMonom> const& coefficients,
                                    In const& monomialValues,
                                    Out& polynomialValues)
    {
      for (int i = 0; i < sizePolynom; ++i)
      {
        squeezeTensor(polynomialValues[i]) = 0;
        for (int j = 0; j < sizeMonom; ++j)
          squeezeTensor(polynomialValues[i]) += coefficients[i][j]*monomialValues[j];
      }
    }

    /** \brief Associations of the Hermite degrees of freedom to subentities of the
     * reference simplex
     *
     * \tparam dim Dimension of the reference simplex
     */
    template<int dim, bool reduced>
    class CubicHermiteLocalCoefficients
    {
    public:
      using size_type = std::size_t;

      CubicHermiteLocalCoefficients()
        : localKeys_(size())
      {
        static_assert((dim > 0) and (dim <= 3), "CubicHermiteLocalCoefficients only implemented for dim=1,2,3");
        static_assert((not reduced) or (dim == 2), "Reduced version of CubicHermiteLocalCoefficients only implemented for dim=2");
        for (size_type i = 0; i < (dim +1); ++i)
        {
          // dim derivatives + 1 evaluation dofs per vertex
          for (size_type k = 0; k < (dim +1); ++k)
            localKeys_[(dim +1) * i + k] = LocalKey(i, dim, k);
        }
        if constexpr (not reduced)
        {
          // 1 evaluation per element (2d) / facets (3d)
          for (size_type i = 0; i < (dim - 1) * (dim - 1); ++i)
            localKeys_[(dim +1) * (dim +1) + i] = LocalKey(i, (dim == 2) ? 0 : 1, 0); // inner dofs
        }
      }

      /** \brief number of coefficients
       */
      static constexpr size_type size()
      {
        if constexpr (dim==1)
          return 4;
        if constexpr ((dim==2) and (reduced))
          return 9;
        if constexpr ((dim==2) and (not reduced))
          return 10;
        if constexpr (dim==3)
          return 20;
        return 0;
      }

      /** \brief get i'th index
       */
      LocalKey const &localKey(size_type i) const
      {
        return localKeys_[i];
      }

    private:
      std::vector<LocalKey> localKeys_;
    };


    /** \brief Implementation of hermite Polynomials
     * \tparam D Type to represent the field in the domain
     * \tparam R Type to represent the field in the range
     * \tparam dim Dimension of the domain simplex
     */
    template<class D, class R, int dim, bool reduced>
    class CubicHermiteReferenceLocalBasis
    {
    public:
      using Traits = H2LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>,
                                        FieldMatrix<R, 1, dim>, FieldMatrix<R, dim, dim>>;

    private:

      /**
       * \brief Get the Hermite Coefficients Matrix
       * \return FieldMatrix<F, (possibly reduced) size, size>
       *  where size is the dimension of the cubic polynomial space
       *
       * This returns the basis transformation matrix from a monomial
       * basis to the Hermite basis on the reference domain.
       * I.e. the i-th row of the returned matrix contains the
       * coefficients of the i-th Hermite basis function with respect to a
       * monomial basis. The monomials are enumerated as in the
       * MonomialSet of the corresponding dimension and order.
       */
      static constexpr auto getCubicHermiteCoefficients()
      {
        if constexpr (dim == 1)
          return Dune::FieldMatrix<D, 4, 4>({{1, 0, -3, 2}, {0, 1, -2, 1}, {0, 0, 3, -2}, {0, 0, -1, 1}});
        else if constexpr (dim == 2)
        {
          if constexpr (reduced) {
            auto w = std::array<D, 9>{1. / 3,  1. / 18, 1. / 18, 1. / 3, -1. / 9,
                                      1. / 18, 1. / 3,  1. / 18, -1. / 9};
            return Dune::FieldMatrix<D, 9, 10>({
                {1, 0, 0, -3, -13 + w[0] * 27, -3,  2, 13 - w[0] * 27, 13 - w[0] * 27,  2},
                {0, 1, 0, -2,  -3 + w[1] * 27,  0,  1,  3 - w[1] * 27,  2 - w[1] * 27,  0},
                {0, 0, 1,  0,  -3 + w[2] * 27, -2,  0,  2 - w[2] * 27,  3 - w[2] * 27,  1},
                {0, 0, 0,  3,  -7 + w[3] * 27,  0, -2,  7 - w[3] * 27,  7 - w[3] * 27,  0},
                {0, 0, 0, -1,   2 + w[4] * 27,  0,  1, -2 - w[4] * 27, -2 - w[4] * 27,  0},
                {0, 0, 0,  0,  -1 + w[5] * 27,  0,  0,  2 - w[5] * 27,  1 - w[5] * 27,  0},
                {0, 0, 0,  0,  -7 + w[6] * 27,  3,  0,  7 - w[6] * 27,  7 - w[6] * 27, -2},
                {0, 0, 0,  0,  -1 + w[7] * 27,  0,  0,  1 - w[7] * 27,  2 - w[7] * 27,  0},
                {0, 0, 0,  0,   2 + w[8] * 27, -1,  0, -2 - w[8] * 27, -2 - w[8] * 27,  1},
            });
          }
          else
            return Dune::FieldMatrix<D, 10,10>({
                                        {1, 0, 0, -3, -13, -3,  2,  13,  13,  2},
                                        {0, 1, 0, -2,  -3,  0,  1,   3,   2,  0},
                                        {0, 0, 1,  0,  -3, -2,  0,   2,   3,  1}, // l_2
                                        {0, 0, 0,  3,  -7,  0, -2,   7,   7,  0},
                                        {0, 0, 0, -1,   2,  0,  1,  -2,  -2,  0},
                                        {0, 0, 0,  0,  -1,  0,  0,   2,   1,  0},
                                        {0, 0, 0,  0,  -7,  3,  0,   7,   7, -2}, // l_6
                                        {0, 0, 0,  0,  -1,  0,  0,   1,   2,  0},
                                        {0, 0, 0,  0,   2, -1,  0,  -2,  -2,  1},
                                        {0, 0, 0,  0,  27,  0,  0, -27, -27,  0}}); // l_9, inner dof
        }
        else if constexpr (dim == 3)
        {
          return Dune::FieldMatrix<D, 20,20>({{1, 0,  0,  0, -3, -13, -3, -13, -13, -3, // deg 0 to 2
                                      2, 13, 13, 2, 13, 33,  13, 13,  13,  2}, // deg 3
                                      {0, 1, 0, 0,/*xx*/ -2, /*xy*/-3,/*yy*/ 0,/*xz*/ -3,/*yz*/ 0,/*zz*/ 0, 1, 3, 2, 0, 3, 4, 0, 2, 0, 0},
                                      {0, 0, 1, 0, 0, -3, -2, 0, -3, 0, 0, 2, 3, 1, 0, 4, 3, 0, 2, 0},
                                      {0, 0, 0, 1, 0, 0, 0, -3, -3, -2, 0, 0, 0, 0, 2, 4, 2, 3, 3, 1},
                                      {0,  0, 0, 0, 3, -7, 0, -7, 0, 0, // l_4
                                      -2, 7, 7, 0, 7, 7,  0, 7,  0, 0},
                                      {0, 0, 0, 0, -1, 2, 0, 2, 0, 0, 1, -2, -2, 0, -2, -2, 0, -2, 0, 0},
                                      {0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0},
                                      {0, 0, 0, 0,  0, -7, 3, 0, -7, 0, // l_8
                                      0, 7, 7, -2, 0, 7,  7, 0, 7,  0},
                                      {0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 2, -1, 0, 2, 0, 0, -2, -2, 1, 0, -2, -2, 0, -2, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0},
                                      {0, 0, 0, 0, 0, 0, 0, -7, -7, 3, // l_12
                                      0, 0, 0, 0, 7, 7, 7, 7,  7,  -2},
                                      {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 2, 2, -1, 0, 0, 0, 0, -2, -2, -2, -2, -2, 1},
                                      // l_16, from here on inner dofs
                                      {0, 0,   0,   0, 0, 27,  0, 0, 0, 0, // bottom
                                      0, -27, -27, 0, 0, -27, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0,   0,   0, 27,  0, 0, // front
                                      0, 0, 0, 0, -27, -27, 0, -27, 0, 0},
                                      {0, 0, 0, 0, 0, 0,   0,   0, 27,  0, // left
                                      0, 0, 0, 0, 0, -27, -27, 0, -27, 0},
                                      {0, 0, 0, 0, 0, 0,  0, 0, 0, 0, // right
                                      0, 0, 0, 0, 0, 27, 0, 0, 0, 0}});
        }
      }

      static constexpr auto referenceBasisCoefficients = getCubicHermiteCoefficients();
      static constexpr MonomialSet<typename Traits::RangeFieldType, dim, 3> monomials = {};

    public:

      CubicHermiteReferenceLocalBasis()
      {
        static_assert((dim > 0) and (dim <= 3), "CubicHermiteReferenceLocalBasis only implemented for dim=1,2,3");
        static_assert((not reduced) or (dim == 2), "Reduced version of CubicHermiteReferenceLocalBasis only implemented for dim=2");
      }

      /** The number of basis functions in the basis
       */
      static constexpr unsigned int size()
      {
        return CubicHermiteLocalCoefficients<dim,reduced>::size();
      }

      /** The polynomial order of the basis
       */
      unsigned int order() const
      {
        return 3;
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
    template<class D, int dim, bool reduced = false>
    class CubicHermiteLocalInterpolation
    {
      using size_type = std::size_t;

      static constexpr unsigned int size()
      {
        return CubicHermiteLocalCoefficients<dim,reduced>::size();
      }

      using FunctionalDescriptor = Dune::Functions::Impl::FunctionalDescriptor<dim>;

    public:

      CubicHermiteLocalInterpolation()
      {
        static_assert((dim > 0) and (dim <= 3), "CubicHermiteLocalInterpolation only implemented for dim=1,2,3");
        static_assert((not reduced) or (dim == 2), "Reduced version of CubicHermiteLocalInterpolation only implemented for dim=2");
        if constexpr (dim==1)
        {
          descriptors_[0] = FunctionalDescriptor();
          descriptors_[1] = FunctionalDescriptor({1});
          descriptors_[2] = FunctionalDescriptor();
          descriptors_[3] = FunctionalDescriptor({1});
        }
        if constexpr (dim==2)
        {
          descriptors_[0] = FunctionalDescriptor();
          descriptors_[1] = FunctionalDescriptor({1,0});
          descriptors_[2] = FunctionalDescriptor({0,1});
          descriptors_[3] = FunctionalDescriptor();
          descriptors_[4] = FunctionalDescriptor({1,0});
          descriptors_[5] = FunctionalDescriptor({0,1});
          descriptors_[6] = FunctionalDescriptor();
          descriptors_[7] = FunctionalDescriptor({1,0});
          descriptors_[8] = FunctionalDescriptor({0,1});
          if (not reduced)
            descriptors_[9] = FunctionalDescriptor();
        }
        if constexpr (dim==3)
        {
          descriptors_[0]  = FunctionalDescriptor();
          descriptors_[1]  = FunctionalDescriptor({1,0,0});
          descriptors_[2]  = FunctionalDescriptor({0,1,0});
          descriptors_[3]  = FunctionalDescriptor({0,0,1});
          descriptors_[4]  = FunctionalDescriptor();
          descriptors_[5]  = FunctionalDescriptor({1,0,0});
          descriptors_[6]  = FunctionalDescriptor({0,1,0});
          descriptors_[7]  = FunctionalDescriptor({0,0,1});
          descriptors_[8]  = FunctionalDescriptor();
          descriptors_[9]  = FunctionalDescriptor({1,0,0});
          descriptors_[10] = FunctionalDescriptor({0,1,0});
          descriptors_[11] = FunctionalDescriptor({0,0,1});
          descriptors_[12] = FunctionalDescriptor();
          descriptors_[13] = FunctionalDescriptor({1,0,0});
          descriptors_[14] = FunctionalDescriptor({0,1,0});
          descriptors_[15] = FunctionalDescriptor({0,0,1});
          descriptors_[16] = FunctionalDescriptor();
          descriptors_[17] = FunctionalDescriptor();
          descriptors_[18] = FunctionalDescriptor();
          descriptors_[19] = FunctionalDescriptor();
        }
      }

      /** \brief bind the Interpolation to an element and a local state.
       */
      template<class Element>
      void bind( Element const &element, std::array<D, dim+1>const& averageVertexMeshSize)
      {
        averageVertexMeshSize_ = &averageVertexMeshSize;
      }

      /** \brief Evaluate a given function and its derivatives at the nodes
       *
       * \tparam F Type of function to evaluate
       * \tparam C Type used for the values of the function
       * \param[in] f Function to evaluate
       * \param[out] out Array of function values
       */
      template<class F, class C>
      void interpolate(const F &f, std::vector<C> &out) const
      {
        out.resize(size());
        auto df = derivative(f);
        auto const &refElement = Dune::ReferenceElements<D, dim>::simplex();

        // Iterate over vertices, dim derivative +1 evaluation dofs per vertex
        for (int i = 0; i < (dim+1); ++i)
        {
          auto x = refElement.position(i, dim);
          auto&& derivativeValue = df(x);
          out[i * (dim +1)] = f(x);
          for (int d = 0; d < dim; ++d)
            out[i * (dim+1) + d + 1] = squeezeTensor(derivativeValue)[d] * (*averageVertexMeshSize_)[i];
        }

        if constexpr (not reduced)
        {
          for (size_type i = 0; i < (dim - 1) * (dim - 1); ++i)
            out[(dim +1) * (dim +1) + i] = f(refElement.position(i, (dim == 2) ? 0 : 1));
        }
      }

      /**
       * \brief Get the object that describes which type of evaluations is performed by the dual basis
       */
      const FunctionalDescriptor& functionalDescriptor(size_type i) const
      {
        return descriptors_[i];
      }

    protected:
      std::array<D, dim+1> const* averageVertexMeshSize_;
      std::array<FunctionalDescriptor, size()> descriptors_;
    };

    template<class D, class R, int dim , bool reduced>
    struct CubicHermiteLocalBasisTraits
      : public H2LocalBasisTraits<D, dim, Dune::FieldVector<D,dim>, R, 1,
              Dune::FieldVector<R,1>, Dune::FieldMatrix<R,1,dim>, Dune::FieldMatrix<R,dim,dim>>
    {};

    /** \brief Hermite finite element for simplices, as defined on the reference Element.
     * Note, that this is a non affine-equivalent finite element, that requires an additional transformation to the relate reference basis with the pullbacks of global basis.
     * For more Details, see <dune/functions/functionspacebases/hermitebasis.hh>.
     *
     * \tparam D Type used for domain coordinates
     * \tparam R Type used for function values
     * \tparam dim dimension of the reference element
     */
    template<class D, class R, int dim, bool reduced = false>
    class CubicHermiteLocalFiniteElement
      : public Impl::TransformedFiniteElementMixin<CubicHermiteLocalFiniteElement<D,R,dim,reduced>, CubicHermiteLocalBasisTraits<D, R, dim, reduced>>
    {
      using Base = Impl::TransformedFiniteElementMixin< CubicHermiteLocalFiniteElement<D,R,dim,reduced>, CubicHermiteLocalBasisTraits<D, R, dim, reduced>>;
      friend class Impl::TransformedLocalBasis<CubicHermiteLocalFiniteElement<D,R,dim,reduced>, CubicHermiteLocalBasisTraits<D, R, dim, reduced>>;

    public:

      CubicHermiteLocalFiniteElement()
        : Base()
      {
        static_assert((dim > 0) and (dim <= 3), "CubicHermiteLocalFiniteElement only implemented for dim=1,2,3");
        static_assert((not reduced) or (dim == 2), "Reduced version of CubicHermiteLocalFiniteElement only implemented for dim=2");
      }

      /** \brief Export number types, dimensions, etc.
       */
      using size_type = std::size_t;
      using Traits = LocalFiniteElementTraits<
          Impl::TransformedLocalBasis<CubicHermiteLocalFiniteElement<D,R,dim,reduced>, CubicHermiteLocalBasisTraits<D, R, dim, reduced>>,
          Impl::CubicHermiteLocalCoefficients<dim, reduced>,
          Impl::CubicHermiteLocalInterpolation<D, dim, reduced>>;

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
        return Impl::CubicHermiteLocalCoefficients<dim,reduced>::size();
      }

      /** Binds the Finite Element to an element.
       */
      template<class Mapper, class Element>
      void bind(Mapper const& vertexMapper, std::vector<D> const& globalAverageVertexMeshSize, Element const &e)
      {
        // Cache average mesh size for each vertex
        for (auto i : range(dim+1))
          averageVertexMeshSize_[i] = globalAverageVertexMeshSize[vertexMapper.subIndex(e, i, dim)];

        // Bind LocalInterpolation to updated local state
        interpolation_.bind(e, averageVertexMeshSize_);

        // Compute local transformation matrices for each vertex
        const auto& geometry = e.geometry();
        const auto& refElement = Dune::ReferenceElements<typename Element::Geometry::ctype, dim>::simplex();
        for (auto i : range(dim+1))
        {
          scaledVertexJacobians_[i] = geometry.jacobian(refElement.position(i, dim));
          scaledVertexJacobians_[i] /= averageVertexMeshSize_[i];
        }
      }

    protected:

      /** \brief Returns the local basis, i.e., the set of shape functions
       */
      Impl::CubicHermiteReferenceLocalBasis<D, R, dim, reduced> const& referenceLocalBasis() const
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
        assert(inValues.size() == size());
        assert(outValues.size() == inValues.size());
        auto inIt = inValues.begin();
        auto outIt = outValues.begin();

        for (auto vertex : Dune::range((dim +1)))
        {
          *outIt = *inIt; // value dof is not transformed
          outIt++, inIt++;
          // transform the gradient dofs together
          for (auto &&[row_i, i] : sparseRange(scaledVertexJacobians_[vertex]))
          {
            outIt[i] = 0.;
            for (auto &&[val_i_j, j] : sparseRange(row_i))
              outIt[i] += val_i_j * inIt[j];
          }
          // increase pointer by size of gradient = dim
          outIt += dim, inIt += dim;
        }

        // For the non-reduced case: Copy all remaining inner dofs
        if constexpr (dim > 1 and (not reduced))
          std::copy(inIt, inValues.end(), outIt);
      }

    private:

      typename Impl::CubicHermiteReferenceLocalBasis<D, R, dim, reduced> basis_;
      typename Traits::LocalCoefficientsType coefficients_;
      typename Traits::LocalInterpolationType interpolation_;
      // the transformation to correct the lack of affine equivalence boils down to
      // one transformation matrix per vertex
      std::array<Dune::FieldMatrix<R, dim, dim>, dim+1> scaledVertexJacobians_;
      // the local state, i.e. a collection of global information restricted to this element
      std::array<D, dim+1> averageVertexMeshSize_;

    };

  } // end namespace Impl



  // *****************************************************************************
  // This is the reusable part of the basis. It contains
  //
  //   CubicHermitePreBasis
  //   CubicHermiteNode
  //
  // The pre-basis allows to create the others and is the owner of possible shared
  // state. These components do _not_ depend on the global basis and local view
  // and can be used without a global basis.
  // *****************************************************************************

  template<class GV, class R, bool reduced>
  class CubicHermiteNode
    : public LeafBasisNode
  {
    using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;
    using FiniteElement = typename Impl::CubicHermiteLocalFiniteElement<typename GV::ctype, R, GV::dimension, reduced>;

    CubicHermiteNode(Mapper const& m, std::vector<typename GV::ctype> const& averageVertexMeshSize)
      : element_(nullptr)
      , vertexMapper_(&m)
      , averageVertexMeshSize_(&averageVertexMeshSize)
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
      finiteElement_.bind(*vertexMapper_, *averageVertexMeshSize_, *element_);
      this->setSize(finiteElement_.size());
    }

    //! The order of the local basis.
    unsigned int order() const { return finiteElement_.localBasis().order(); }

  protected:
    FiniteElement finiteElement_;
    Element const* element_;
    Mapper const* vertexMapper_;
    std::vector<typename GV::ctype> const* averageVertexMeshSize_;
  };


  /**
   * \brief A pre-basis for a Hermitebasis
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV  The grid view that the FE basis is defined on
   * \tparam R   Range type used for shape function values
   * \note This only works for simplex grids
   */
  template<class GV, class R, bool reduced = false>
  class CubicHermitePreBasis
    : public LeafPreBasisMapperMixin<GV>
  {
    using Base = LeafPreBasisMapperMixin<GV>;
    using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
    using Element = typename GV::template Codim<0>::Entity;
    using D = typename GV::ctype;
    static const std::size_t dim = GV::dimension;

    // helper methods to assign each subentity the number of dofs. Used by the LeafPreBasisMapperMixin.
    static constexpr auto cubicHermiteMapperLayout(Dune::GeometryType type, int gridDim)
    {
      if (type.isVertex())
        return 1 + gridDim; // one evaluation dof and gridDim derivative dofs per vertex
      if (gridDim == 1)    // in 1d there are no other dofs
        return 0;
      // in 2d we have one inner dof (i.e. on the triangle) or non for the reduced case
      // and in 3d we have one dof on each face (i.e. on each triangle)
      if ((type.isTriangle()) and (not reduced))
        return 1;
      else
        return 0; // this case is only entered for the interior of the 3d element. There are no dofs.
    }

  public:
    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

    //! Template mapping root tree path to type of created tree node
    using Node = CubicHermiteNode<GridView, R,reduced>;

  public:

    //! Constructor for a given grid view object
    CubicHermitePreBasis(const GV &gv)
      : Base(gv, cubicHermiteMapperLayout)
      , vertexMapper_({gv, mcmgVertexLayout()})
    {
      static_assert((dim > 0) and (dim <= 3), "CubicHermitePreBasis only implemented for dim=1,2,3");
      static_assert((not reduced) or (dim == 2), "Reduced version of CubicHermitePreBasis only implemented for dim=2");
      averageVertexMeshSize_ = Impl::computeAverageSubEntityMeshSize<D>(vertexMapper_);
    }

    //! Update the stored grid view, to be called if the grid has changed
    void update(GridView const &gv)
    {
      Base::update(gv);
      vertexMapper_.update(this->gridView());
      averageVertexMeshSize_ = Impl::computeAverageSubEntityMeshSize<D>(vertexMapper_);
    }

    /**
     * \brief Create tree node
     */
    Node makeNode() const
    {
      return Node{vertexMapper_, averageVertexMeshSize_};
    }

  protected:

    SubEntityMapper vertexMapper_;
    std::vector<D> averageVertexMeshSize_;

  }; // class CubicHermitePreBasis

  namespace BasisFactory
  {

    /**
     * \brief construct a PreBasisFactory for the full cubic Hermite Finite Element
     *
     * \tparam R RangeFieldType
     * \return the PreBasisFactory
     * \relates CubicHermiteBasis
     */
    template<class R = double>
    auto cubicHermite()
    {
      return [=](auto const &gridView) {
        return CubicHermitePreBasis<std::decay_t<decltype(gridView)>, R>(gridView);
      };
    }

    /**
     * \brief construct a PreBasisFactory for the reduced cubic Hermite Finite Element
     *
     * \tparam R RangeFieldType
     * \return the PreBasisFactory
     * \relates CubicHermiteBasis
     */
    template<class R = double>
    auto reducedCubicHermite()
    {
      return [=](auto const &gridView) {
        return CubicHermitePreBasis<std::decay_t<decltype(gridView)>, R, true>(gridView);
      };
    }

  } // end namespace BasisFactory


} // end namespace Dune::Functions

#endif
