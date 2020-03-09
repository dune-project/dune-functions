// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GLOBALVALUEDLOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GLOBALVALUEDLOCALFINITEELEMENT_HH

#include <array>
#include <numeric>

#include <dune/common/deprecated.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune::Functions::Impl
{

  /** \brief Transforms shape function values and derivatives from reference element coordinates
   *   to world coordinates using the contravariant Piola transform
   */
  struct ContravariantPiolaTransformator
  {
    template<typename Values, typename LocalCoordinate, typename Geometry>
    auto apply(Values&& valuesLocal,
      const LocalCoordinate& xi,
      const Geometry& geometry)
    {
      auto jacobianTransposed = geometry.jacobianTransposed(xi);
      auto integrationElement = geometry.integrationElement(xi);

      for (auto& value : valuesLocal)
      {
        auto tmp = value;
        jacobianTransposed.mtv(tmp, value);
        value /= integrationElement;
      }

      return std::move(valuesLocal);
    }
  };


  template<class LocalValuedLocalBasis, class Element>
  class GlobalValuedLocalBasis
  {
  public:
    using Traits = typename LocalValuedLocalBasis::Traits;

    void bind(const LocalValuedLocalBasis* localValuedLocalBasis, const Element* element)
    {
      localValuedLocalBasis_ = localValuedLocalBasis;
      element_ = element;
    }

    /** \brief Number of shape functions
     */
    auto size() const
    {
      return localValuedLocalBasis_->size();
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      localValuedLocalBasis_->evaluateFunction(x,out);

      // Apply the Piola transform
      auto jacobianTransposed = element_->geometry().jacobianTransposed(x);
      auto integrationElement = element_->geometry().integrationElement(x);

      for (auto& value : out)
      {
        auto tmp = value;
        jacobianTransposed.mtv(tmp, value);
        value /= integrationElement;
      }
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference element where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      localValuedLocalBasis_->evaluateJacobian(x,out);
      DUNE_THROW(NotImplemented, "Piola transform is missing in 'evaluateJacobian'!");
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int,2>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      localValuedLocalBasis_->partial(order, in, out);
      DUNE_THROW(NotImplemented, "Piola transform is missing in 'partial'!");
    }

    //! \brief Polynomial order of the shape functions
    auto order() const
    {
      return localValuedLocalBasis_->order();
    }

    const LocalValuedLocalBasis* localValuedLocalBasis_;
    const Element* element_;
  };

  /** \brief Lagrange finite element for simplices with arbitrary compile-time dimension and polynomial order
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   */
  template<class LocalValuedLFE, class Element>
  class GlobalValuedLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<GlobalValuedLocalBasis<typename LocalValuedLFE::Traits::LocalBasisType,Element>,
                                            typename LocalValuedLFE::Traits::LocalCoefficientsType,
                                            typename LocalValuedLFE::Traits::LocalInterpolationType>;

    GlobalValuedLocalFiniteElement() {}

    void bind(const LocalValuedLFE& localValuedLFE, const Element* element)
    {
      globalValuedLocalBasis_.bind(&localValuedLFE.localBasis(), element);
      localValuedLFE_ = &localValuedLFE;
    }

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis() const
    {
      return globalValuedLocalBasis_;
    }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
      return localValuedLFE_->localCoefficients();
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType& localInterpolation() const
    {
      // TODO: Modify this!
      return localValuedLFE_->localInterpolation();
    }

    /** \brief The number of shape functions */
    std::size_t size() const
    {
      return localValuedLFE_->size();
    }

    /** \brief The reference element that the local finite element is defined on
     */
    GeometryType type() const
    {
      return localValuedLFE_->type();
    }

  private:

    GlobalValuedLocalBasis<typename LocalValuedLFE::Traits::LocalBasisType,Element> globalValuedLocalBasis_;
    const LocalValuedLFE* localValuedLFE_;
  };

}        // namespace Dune::Functions::Impl

#endif   // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GLOBALVALUEDLOCALFINITEELEMENT_HH
