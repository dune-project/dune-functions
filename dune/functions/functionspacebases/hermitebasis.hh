// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HERMITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HERMITEBASIS_HH

#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <type_traits>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/transformedfiniteelementmixin.hh>
#include <dune/functions/analyticfunctions/monomialset.hh>


namespace Dune
{
namespace Functions
{
  template<class DF, int n, class D, class RF, int m, class R, class J, class H>
  struct H2LocalBasisTraits : public LocalBasisTraits<DF, n, D, RF, m, R, J> {
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
  /**
   * \brief Multiply the evaluations of the monomialSet (see dune/functions/analyticfunctions/monomialset.hh) with a coefficient matrix, here FieldMatrix.
   * \tparam KCoeff The Field type of the coefficient matrix.
   * \tparam KMonom The range type of the monomials.
   * \tparam sizePolynom The number of polynomials to evaluate
   * \tparam sizeMonom The number of monomials handed as input.
   * \note This function also turns the return types from the dune-functions::DifferentiableFunction interface into those of the dune-localfunctions::LocalBasis interface.
   */
  template<class KCoeff, int sizePolynom, class KMonom, int sizeMonom>
  void multiplyWithCoefficentMatrix(Dune::FieldMatrix<KCoeff, sizePolynom, sizeMonom> const& coefficients,
                                   Dune::FieldVector<KMonom, sizeMonom> const& monomialValues,
                                   std::vector<Dune::FieldVector<typename Dune::PromotionTraits<KCoeff, KMonom>::PromotedType, 1>>& polynomialValues)
  {

    polynomialValues.resize(sizePolynom);
    std::fill(std::begin(polynomialValues), std::end(polynomialValues), 0.);
    for (auto&& i : Dune::range(sizePolynom))
      for (auto&& j : Dune::range(sizeMonom))
        polynomialValues[i][0] += coefficients[i][j]*monomialValues[j];
  }

  template<class KCoeff, int sizePolynom, class KMonom, int sizeMonom, int dim>
  void multiplyWithCoefficentMatrix(Dune::FieldMatrix<KCoeff, sizePolynom, sizeMonom> const& coefficients,
                                   Dune::FieldMatrix<KMonom, sizeMonom, dim> const& monomialJacobians,
                                   std::vector<Dune::FieldMatrix<typename Dune::PromotionTraits<KCoeff, KMonom>::PromotedType, 1,dim>>& polynomialJacobians)
  {
    polynomialJacobians.resize(sizePolynom);
    std::fill(std::begin(polynomialJacobians), std::end(polynomialJacobians), 0.);
    for (auto&& i : Dune::range(sizePolynom))
      for (auto&& j : Dune::range(sizeMonom))
        polynomialJacobians[i][0] += coefficients[i][j]*monomialJacobians[j];
  }

  template<class KCoeff, int sizePolynom, class KMonom, int sizeMonom, std::size_t sizeMonom2, int dim>
  void multiplyWithCoefficentMatrix(Dune::FieldMatrix<KCoeff, sizePolynom, sizeMonom> const& coefficients,
                                   std::array<Dune::FieldMatrix<KMonom, dim, dim>, sizeMonom2> const& monomialHessians,
                                   std::vector<Dune::FieldMatrix<typename Dune::PromotionTraits<KCoeff, KMonom>::PromotedType, dim, dim>>& polynomialHessians)
  {
    polynomialHessians.resize(sizePolynom);
    std::fill(std::begin(polynomialHessians), std::end(polynomialHessians), 0.);
    for (auto&& i : Dune::range(sizePolynom))
      for (auto&& j : Dune::range(sizeMonom))
        polynomialHessians[i] += coefficients[i][j]*monomialHessians[j];
  }

  /**
   * \brief Implementation of hermite Polynomials
   * \tparam D Type to represent the field in the domain
   * \tparam R Type to represent the field in the range
   * \tparam dim Dimension of the domain simplex
   */
  template<class D, class R, int dim, bool reduced>
  class HermiteLocalBasis
  {
    public:
      using Traits = H2LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>,
                                        FieldMatrix<R, 1, dim>, FieldMatrix<R, dim, dim>>;

    private:
      /**
       * @brief Get the Hermite Coefficients Matrix
       * @return FieldMatrix<F, (possibly reduced) size, size>
       *  where size is the dimension of the cubic polynomial space
       */
      static constexpr auto getHermiteCoefficients()
      {
        static_assert(dim > 0 and dim < 4 and not(reduced and dim != 2));

        if constexpr (dim == 1) {
          return Dune::FieldMatrix<D, 4, 4>({{1, 0, -3, 2}, {0, 1, -2, 1}, {0, 0, 3, -2}, {0, 0, -1, 1}});
        } else if constexpr (dim == 2) {
          if constexpr (reduced) {
            auto w = std::array<D, 9>{1. / 3,  1. / 18, 1. / 18, 1. / 3, -1. / 9,
                                      1. / 18, 1. / 3,  1. / 18, -1. / 9};
            return Dune::FieldMatrix<D, 9, 10>({
                {1, 0, 0, -3, -13 + w[0] * 27, -3, 2, 13 - w[0] * 27, 13 - w[0] * 27, 2},
                {0, 1, 0, -2, -3 + w[1] * 27, 0, 1, 3 - w[1] * 27, 2 - w[1] * 27, 0},
                {0, 0, 1, 0, -3 + w[2] * 27, -2, 0, 2 - w[2] * 27, 3 - w[2] * 27, 1},
                {0, 0, 0, 3, -7 + w[3] * 27, 0, -2, 7 - w[3] * 27, 7 - w[3] * 27, 0},
                {0, 0, 0, -1, 2 + w[4] * 27, 0, 1, -2 - w[4] * 27, -2 - w[4] * 27, 0},
                {0, 0, 0, 0, -1 + w[5] * 27, 0, 0, 2 - w[5] * 27, 1 - w[5] * 27, 0},
                {0, 0, 0, 0, -7 + w[6] * 27, 3, 0, 7 - w[6] * 27, 7 - w[6] * 27, -2},
                {0, 0, 0, 0, -1 + w[7] * 27, 0, 0, 1 - w[7] * 27, 2 - w[7] * 27, 0},
                {0, 0, 0, 0, 2 + w[8] * 27, -1, 0, -2 - w[8] * 27, -2 - w[8] * 27, 1},
            });
          } else
            return Dune::FieldMatrix<D, 10,10>({{1, 0, 0, -3, -13, -3, 2, 13, 13, 2},
                                        {0, 1, 0, -2, -3, 0, 1, 3, 2, 0},
                                        {0, 0, 1, 0, -3, -2, 0, 2, 3, 1}, // l_2
                                        {0, 0, 0, 3, -7, 0, -2, 7, 7, 0},
                                        {0, 0, 0, -1, 2, 0, 1, -2, -2, 0},
                                        {0, 0, 0, 0, -1, 0, 0, 2, 1, 0},
                                        {0, 0, 0, 0, -7, 3, 0, 7, 7, -2}, // l_6
                                        {0, 0, 0, 0, -1, 0, 0, 1, 2, 0},
                                        {0, 0, 0, 0, 2, -1, 0, -2, -2, 1},
                                        {0, 0, 0, 0, 27, 0, 0, -27, -27, 0}}); // l_9, inner dof
        } else if constexpr (dim == 3) {
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

      static constexpr auto referenceBasisCoefficients = getHermiteCoefficients();
      MonomialSet<typename Traits::RangeFieldType, dim, 3> monomials;

    public:
      static_assert(not reduced || dim == 2, "Reduced Hermite element only implemented for 2d");
      static constexpr int coeffSize = (dim == 1)   ? 4
                                                : (dim == 2) ? ((reduced) ? 9 : 10)
                                                            : 20;
      HermiteLocalBasis()
      {
        if (not (dim <= 3))
          DUNE_THROW(Dune::NotImplemented, "only implemented for dim <= 3");
      }

      /** The number of basis functions in the basis
       */
      static constexpr unsigned int size() { return coeffSize; }

      /** The polynomial order of the basis
       */
      unsigned int order() const { return 3; }

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
        else if (totalOrder == 1){
          evaluateJacobian(in,jacobiansBuffer_);
          std::size_t which = std::max_element(order.begin(), order.end()) - order.begin();
          for (auto i : Dune::range(size()))
            out[i] = jacobiansBuffer_[i][0][which];
        }
        else if (totalOrder == 2){
          evaluateHessian(in, hessianBuffer_);
          std::size_t first, second;
          first = std::max_element(order.begin(), order.end()) - order.begin();
          if (order[first] == 2){
            second = first;
          } else {
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

  /** \brief Associations of the Hermite degrees of freedom to subentities of the
   * reference simplex
   *
   * \tparam dim Dimension of the reference simplex
   */
  template<unsigned int dim, bool reduced>
  class HermiteLocalCoefficients
  {
    public:
      using size_type = std::size_t;

    private:
      static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1;

    public:
      HermiteLocalCoefficients() : localKeys_(size())
      {
        static_assert(dim <= 3, "HermiteLocalCoefficients only implemented for dim<=3!");
        size_type numberOfVertices = dim + 1;
        size_type numberOfInnerDofs = (dim - 1) * (dim - 1); // probably incorrect for dim > 3
        for (size_type i = 0; i < numberOfVertices; ++i)     // subentities: vertices
        {
          for (size_type k = 0; k < numberOfVertices; ++k) // dim + 1 dofs per subentity
            localKeys_[numberOfVertices * i + k] = LocalKey(i, dim, k);
        }
        if constexpr (not reduced)
          for (size_type i = 0; i < numberOfInnerDofs; ++i) // subentities: element
            localKeys_[numberOfVertices * numberOfVertices + i] =
                LocalKey(i, innerDofCodim, 0); // inner dofs
      }

      /** number of coefficients
       */
      static constexpr size_type size()
      {
        return dim == 1 ? 4 : dim == 2 ? ((reduced) ? 9 : 10) : 20;
      }

      /** get i'th index
       */
      LocalKey const &localKey(size_type i) const { return localKeys_[i]; }

    private:
      std::vector<LocalKey> localKeys_;
  };

    /**
     * \brief Class that evaluates the push forwards of the global nodes of a
     * LocalFunction. It stretches the LocalInterpolation interface, because we
     * evaluate the derivatives of f.
     *
     */
    template<class D, int dim, bool reduced = false>
    class HermiteLocalInterpolation
    {
      using size_type = std::size_t;
      static constexpr size_type numberOfVertices = dim + 1;
      static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1; // probably wrong for dim > 3
      static constexpr size_type numberOfInnerDofs =
          (dim - 1) * (dim - 1); // probably wrong for dim > 3
      static constexpr unsigned int coeffSize = (dim == 1)   ? 4
                                                : (dim == 2) ? ((reduced) ? 9 : 10)
                                                            : 20;
    public:
      HermiteLocalInterpolation(){}

      /** \brief bind the Interpolation to an element and a localState.
      */
      template<class Element>
      void bind( Element const &element, std::vector<D>const& localState)
      {
          localState_ = &localState;
      }

      /** \brief Evaluate a given function and its derivatives at the nodes
      *
      * \tparam F Type of function to evaluate
      * \tparam C Type used for the values of the function
      * \param[in] f Function to evaluate
      * \param[out] out Array of function values
      */
      template<typename F, typename C>
      void interpolate(const F &f, std::vector<C> &out) const
      {
          auto df = derivative(f);
          out.resize(coeffSize);

          auto const &refElement = Dune::ReferenceElements<D, dim>::simplex();
          // Iterate over vertices, dim +1 dofs per vertex
          for (int i = 0; i < dim + 1; ++i) {
          auto x = refElement.position(i, dim);

          auto derivativeValue = df(x);
          out[i * numberOfVertices] = f(x);
          for (int d = 0; d < dim; ++d)
              out[i * numberOfVertices + d + 1] = getPartialDerivative(derivativeValue,d) * (*localState_)[i];
          }

          if constexpr (not reduced)
          for (size_type i = 0; i < numberOfInnerDofs; ++i) {
              out[numberOfVertices * numberOfVertices + i] =
                  f(refElement.position(i, innerDofCodim));
          }
      }

    protected:
      template<class DerivativeType, class FieldType>
      FieldType getPartialDerivative(DerivativeType const &df, std::size_t i) const
      {
          DUNE_THROW(Dune::NotImplemented, "Derivative Type is neither FieldMatrix<double,1,d> nor "
                                          "FieldVector<double,d>");
      }

      template<class FieldType, int d>
      FieldType getPartialDerivative(Dune::FieldVector<FieldType, d> const &df, std::size_t i) const
      {
          return df[i];
      }

      template<class FieldType, int d>
      FieldType getPartialDerivative(Dune::FieldMatrix<FieldType, 1, d> const &df,
                                      std::size_t i) const
      {
        if (df.N() == 1)
          return df[0][i];
        else if (df.M() == 1)
          return df[i][0];
        else
          DUNE_THROW(Dune::NotImplemented, "Derivative of scalar function is a matrix!");
      }

      std::vector<D> const* localState_;
  };



  } // namespace Impl

  template<class D, class R, unsigned int dim , bool reduced>
  struct HermiteLocalBasisTraits
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
  template<class D, class R, unsigned int dim, bool reduced = false>
  class HermiteLocalFiniteElement: public Impl::TransformedFiniteElementMixin<HermiteLocalFiniteElement<D,R,dim,reduced>, HermiteLocalBasisTraits<D, R, dim, reduced>>
  {
    using Base = Impl::TransformedFiniteElementMixin< HermiteLocalFiniteElement<D,R,dim,reduced>, HermiteLocalBasisTraits<D, R, dim, reduced>>;
    friend class Impl::TransformedLocalBasis<HermiteLocalFiniteElement<D,R,dim,reduced>, HermiteLocalBasisTraits<D, R, dim, reduced>>;

    static_assert(dim > 0 && dim < 4);
    static_assert(!(reduced && (dim != 2)));
    static constexpr std::size_t numberOfVertices = dim + 1;
    static constexpr std::size_t numberOfInnerDofs = reduced ? 0 : (dim - 1) * (dim - 1);
    static constexpr std::size_t numberOfVertexDofs = numberOfVertices * numberOfVertices;
  public:

    HermiteLocalFiniteElement()
      :Base()
      {}
    /** \brief Export number types, dimensions, etc.
     */
    using LocalState = typename std::vector<D>;
    using size_type = std::size_t;
    using Traits = LocalFiniteElementTraits<
        Impl::TransformedLocalBasis<HermiteLocalFiniteElement<D,R,dim,reduced>, HermiteLocalBasisTraits<D, R, dim, reduced>>,
        Impl::HermiteLocalCoefficients<dim, reduced>,
        Impl::HermiteLocalInterpolation<D, dim, reduced>>;

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
    static constexpr GeometryType type() { return GeometryTypes::simplex(dim); }

    /** The size of the transformed finite element.
     */
    static constexpr size_type size()
    {
      if constexpr (dim == 1)
        return 4;
      else if constexpr (dim == 2) {
        if constexpr (reduced)
          return 9;
        else
          return 10;
      } else // dim == 3
        return 20;
    }

    /** Binds the Finite Element to an element.
      */
    template<class Mapper, class Element>
    void bind(Mapper const& mapper, std::vector<D> const& data, Element const &e)
    {
      for (auto const &index : range(e.subEntities(dim)))
        localState_.push_back(data[mapper.subIndex(e, index, dim)]);

      fillMatrix(e.geometry(), localState_);
      interpolation_.bind(e, localState_);
    }
  protected:
    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    Impl::HermiteLocalBasis<D, R, dim, reduced> const&referenceLocalBasis() const { return basis_; }

    /** Applies the transformation. Note that we do not distinguish for
      * Scalar/Vector/Matrix Type,
      * but only assume the Values to be Elements of a Vectorspace.
      * We assume random access containers. */
    template<typename InputValues, typename OutputValues>
    void transform(InputValues const &inValues, OutputValues &outValues) const
    {
      assert(inValues.size() == numberOfVertexDofs + numberOfInnerDofs);
      assert(reduced || (outValues.size() == inValues.size()));
      auto inIt = inValues.begin();
      auto outIt = outValues.begin();

      for (auto vertex : Dune::range(numberOfVertices)) {
        *outIt = *inIt; // value dof is not transformed
        outIt++, inIt++;
        // transform the gradient dofs together
        for (auto &&[row_i, i] : sparseRange(subMatrices_[vertex])) {
          outIt[i] = 0.;
          for (auto &&[val_i_j, j] : sparseRange(row_i))
            outIt[i] += val_i_j * inIt[j];
        }
        // increase pointer by size of gradient = dim
        outIt += dim, inIt += dim;
      }

      if constexpr (not reduced)
        // copy all remaining inner dofs
        std::copy(inIt, inValues.end(), outIt);
    }

  private:
    /**
      * \brief Fill the transformationmatrix m
      *
      * \tparam Geometry  the Geometry class
      * \param geometry   the geometry of the element we are bound to
      *
      *
      *      |1          0|} repeat for each vertex
      * m =  |  J/h       |}
      *      |0          1|} repeat for each inner dof i.e. (dim-1)^2 times
      *  where h is the mesh size average over the local vertex patch
      */
    template<class Geometry>
    void fillMatrix(Geometry const &geometry, LocalState const &averageSubEntityMeshSize)
    {
      auto const &refElement = Dune::ReferenceElements<typename Geometry::ctype, dim>::simplex();
      for (std::size_t i = 0; i < numberOfVertices; ++i) // dim + 1 vertices
      {
        subMatrices_[i] = geometry.jacobian(refElement.position(i, dim));
        subMatrices_[i] /= averageSubEntityMeshSize[i];
      }
    }

    // a finite element consists of a basis, coeffiecents and an interpolation
    typename Impl::HermiteLocalBasis<D, R, dim, reduced> basis_;
    typename Traits::LocalCoefficientsType coefficients_;
    typename Traits::LocalInterpolationType interpolation_;
    // the transformation to correct the lack of affine equivalence boils down to
    // one transformation matrix per vertex
    std::array<Dune::FieldMatrix<R, dim, dim>, numberOfVertices> subMatrices_;
    // the local state, i.e. a collection of global information restricted to this element
    LocalState localState_;

  };


  namespace Impl
  {

    // Helper function returning an unordered range
    // of global indices associated to the element.
    // This could be implemented cheaper internally in
    // the MCMGMapper by storing a precomputed
    // container of all subsentities addressed by the layout.
    template<class GridView>
    auto subIndexSet(Dune::MultipleCodimMultipleGeomTypeMapper<GridView> const &mapper,
                    const typename GridView::template Codim<0>::Entity &element)
    {
      using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
      using Index = typename Mapper::Index;
      constexpr auto dimension = GridView::dimension;
      auto subIndices = std::vector<Index>();
      auto referenceElement = Dune::referenceElement<double, dimension>(element.type());
      for (auto codim : Dune::range(dimension + 1)) {
        for (auto subEntity : Dune::range(referenceElement.size(codim))) {
          std::size_t c = mapper.layout()(referenceElement.type(subEntity, codim), dimension);
          if (c > 0) {
            std::size_t firstIndex = mapper.subIndex(element, subEntity, codim);
            for (auto j : Dune::range(firstIndex, firstIndex + c)) {
              subIndices.push_back(j);
            }
          }
        }
      }
      return subIndices;
    }

    // Helper function computing an average mesh size per subentity
    // by averaging over the adjacent elements. This only considers
    // the subentities handled by the given mapper and returns a
    // vector of mesh sizes indixed according to the mapper.
    template<class FieldType = double, class Mapper>
    auto computeAverageSubEntityMeshSize(const Mapper& mapper)
    {
      constexpr auto dimension = Mapper::GridView::dimension;

      std::vector<unsigned int> adjacentElements(mapper.size(), 0);
      std::vector<FieldType> subEntityMeshSize(mapper.size(), 0.0);
      for(const auto& element : Dune::elements(mapper.gridView()))
      {
        auto A = element.geometry().volume();
        for(auto i : Impl::subIndexSet(mapper, element))
        {
          subEntityMeshSize[i] += A;
          ++(adjacentElements[i]);
        }
      }
      for(auto i : Dune::range(mapper.size()))
        subEntityMeshSize[i] = std::pow(subEntityMeshSize[i]/adjacentElements[i], 1./dimension);
      return subEntityMeshSize;
    }

  } // namespace Impl

template<typename GV, class R, bool reduced>
class HermiteNode : public LeafBasisNode
{
    using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;

    using FiniteElement = HermiteLocalFiniteElement<typename GV::ctype, R, GV::dimension, reduced>;


    HermiteNode(Mapper const& m, std::vector<typename GV::ctype> const& data)
    : mapper_(&m), data_(&data)
    {
      this->setSize(finiteElement_.size());
    }

    //! Return current element, throw if unbound
    Element const &element() const { return *element_; }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the
     * dune-localfunctions module
     */
    FiniteElement const &finiteElement() const { return finiteElement_; }

    //! Bind to element.
    void bind(Element const &e)
    {
      element_ = &e;
      finiteElement_.bind(*mapper_, *data_, *element_);
    }

    //! The order of the local basis.
    unsigned int order() const { return finiteElement_.localBasis().order(); }

  protected:
    FiniteElement finiteElement_;
    Element const* element_;
    Mapper const* mapper_;
    std::vector<typename GV::ctype> const* data_;
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
  template<typename GV, typename R, bool reduced = false>
  class HermitePreBasis : public LeafPreBasisMapperMixin<GV>
  {
    using Base = LeafPreBasisMapperMixin<GV>;
    using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

    static const std::size_t dim = GV::dimension;
    using Element = typename GV::template Codim<0>::Entity;
    using D = typename GV::ctype;
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
    using Node = HermiteNode<GridView, R,reduced>;

    static constexpr size_type maxMultiIndexSize = 1;
    static constexpr size_type minMultiIndexSize = 1;
    static constexpr size_type multiIndexBufferSize = 1;

  public:
    //! Constructor for a given grid view object
    HermitePreBasis(const GV &gv)
        : Base(gv, cubicHermiteMapperLayout), mapper_({gv, mcmgVertexLayout()})
    {
      updateState(gv);
      if (dim > 3)
        DUNE_THROW(Dune::NotImplemented, "HermitePreBasis only implemented for dim <= 3");
    }

    //! Update the stored grid view, to be called if the grid has changed
    void update(GridView const &gv)
    {
      Base::update(gv);
      updateState(gv);
    }

    /**
    * \brief Create tree node
    */
    Node makeNode() const { return Node{mapper_, data_}; }

  protected:
    void updateState(GridView const &gridView)
    {
      mapper_.update(gridView);
      data_ = Impl::computeAverageSubEntityMeshSize<D>(mapper_);
    }

    SubEntityMapper mapper_;
    std::vector<D> data_;

  }; // class HermitePreBasis

  namespace BasisFactory
  {

    /**
    * \brief construct a PreBasisFactory for the full cubic Hermite Finite Element
    *
    * \tparam R RangeFieldType
    * \return the PreBasisFactory
    */
    template<typename R = double>
    auto hermite()
    {
      return [=](auto const &gridView) {
        return HermitePreBasis<std::decay_t<decltype(gridView)>, R>(gridView);
      };
    }

    /**
    * \brief construct a PreBasisFactory for the reduced cubic Hermite Finite Element
    *
    * \tparam R RangeFieldType
    * \return the PreBasisFactory
    */
    template<typename R = double>
    auto reducedHermite()
    {
      return [=](auto const &gridView) {
        return HermitePreBasis<std::decay_t<decltype(gridView)>, R, true>(gridView);
      };
    }

  } // namespace BasisFactory

} // namespace Functions
} // namespace Dune

#endif
