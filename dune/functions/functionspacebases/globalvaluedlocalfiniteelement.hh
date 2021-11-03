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
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune::Functions::Impl
{

  /** \brief Transforms shape function values and derivatives from reference element coordinates
   *   to world coordinates using the contravariant Piola transform
   *
   * This transformation preserves normal traces of vector fields.
   * It is therefore the canonical transformation for H(div)-conforming finite elements.
   *
   * See for example:
   *   M. Rognes, R. Kirby, A. Logg, "Efficient Assembly of H(div) and H(curl)
   *   conforming finite elements", SIAM J. Sci. Comput., 2009
   *
   * \todo Currently each method computes jacobianTransposed and integrationElement again.
   *   Can we cache this somehow?
   */
  struct ContravariantPiolaTransformator
  {
    /** \brief Piola-transform a set of shape-function values
     *
     * \param[in,out] values The values to be Piola-transformed
     */
    template<typename Values, typename LocalCoordinate, typename Geometry>
    static auto apply(Values& values,
                      const LocalCoordinate& xi,
                      const Geometry& geometry)
    {
      auto jacobianTransposed = geometry.jacobianTransposed(xi);
      auto integrationElement = geometry.integrationElement(xi);

      for (auto& value : values)
      {
        auto tmp = value;
        jacobianTransposed.mtv(tmp, value);
        value /= integrationElement;
      }
    }

    /** \brief Piola-transform a set of shape-function derivatives
     *
     * \param[in,out] gradients The shape function derivatives to be Piola-transformed
     *
     * \bug The current implementation works only for affine geometries.
     *   The Piola transformation for non-affine geometries requires
     *   second derivatives of the geometry, which we don't get
     *   from the dune-grid Geometry interface.
     */
    template<typename Gradients, typename LocalCoordinate, typename Geometry>
    static auto applyJacobian(Gradients& gradients,
                              const LocalCoordinate& xi,
                              const Geometry& geometry)
    {
      auto jacobianTransposed = geometry.jacobianTransposed(xi);
      auto integrationElement = geometry.integrationElement(xi);
      for (auto& gradient : gradients)
      {
        auto tmp = gradient;
        gradient = 0;
        for (size_t k=0; k<gradient.M(); k++)
          for (size_t l=0; l<tmp.N(); l++)
            // Use sparseRange because jacobianTransposed may be a sparse DiagonalMatrix
            for(auto&& [jacobianTransposed_l_j, j] : sparseRange(jacobianTransposed[l]))
              gradient[j][k] += jacobianTransposed_l_j * tmp[l][k];
        gradient /= integrationElement;
      }
    }

    /** \brief Wrapper around a callable that applies the inverse Piola transform
     *
     * The LocalInterpolation implementations in dune-localfunctions expect local-valued
     * functions, but the ones dune-functions expect global-valued ones.  Therefore,
     * we need to stuff the inverse Piola transform between dune-functions and
     * dune-localfunctions, and this is what this class does.
     */
    template<class Function, class LocalCoordinate, class Element>
    class LocalValuedFunction
    {
      const Function& f_;
      const Element& element_;

    public:

      LocalValuedFunction(const Function& f, const Element& element)
      : f_(f), element_(element)
      {}

      auto operator()(const LocalCoordinate& xi) const
      {
        auto&& f = Dune::Impl::makeFunctionWithCallOperator<LocalCoordinate>(f_);
        auto globalValue = f(xi);

        // Apply the inverse Piola transform
        auto jacobianInverseTransposed = element_.geometry().jacobianInverseTransposed(xi);
        auto integrationElement = element_.geometry().integrationElement(xi);

        auto localValue = globalValue;
        jacobianInverseTransposed.mtv(globalValue, localValue);
        localValue *= integrationElement;

        return localValue;
      }
    };
  };

  /** \brief Transforms shape function values and derivatives from reference element coordinates
   *   to world coordinates using the covariant Piola transform
   *
   * This transformation preserves tangential traces of vector fields.
   * It is therefore the canonical transformation for H(curl)-conforming finite elements.
   *
   * See for example:
   *   M. Rognes, R. Kirby, A. Logg, "Efficient Assembly of H(div) and H(curl)
   *   conforming finite elements", SIAM J. Sci. Comput., 2009
   *
   * \todo Currently each method computes jacobianTransposed and integrationElement again.
   *   Can we cache this somehow?
   */
  struct CovariantPiolaTransformator
  {
    /** \brief Piola-transform a set of shape-function values
     *
     * \param[in,out] values The values to be Piola-transformed
     */
    template<typename Values, typename LocalCoordinate, typename Geometry>
    static auto apply(Values& values,
                      const LocalCoordinate& xi,
                      const Geometry& geometry)
    {
      auto jacobianInverseTransposed = geometry.jacobianInverseTransposed(xi);

      for (auto& value : values)
      {
        auto tmp = value;
        jacobianInverseTransposed.mv(tmp, value);
      }
    }

    /** \brief Piola-transform a set of shape-function derivatives
     *
     * \param[in,out] gradients The shape function derivatives to be Piola-transformed
     *
     * \bug The current implementation works only for affine geometries.
     *   The Piola transformation for non-affine geometries requires
     *   second derivatives of the geometry, which we don't get
     *   from the dune-grid Geometry interface.
     */
    template<typename Gradients, typename LocalCoordinate, typename Geometry>
    static auto applyJacobian(Gradients& gradients,
                              const LocalCoordinate& xi,
                              const Geometry& geometry)
    {
      auto jacobianInverseTransposed = geometry.jacobianInverseTransposed(xi);

      for (auto& gradient : gradients)
      {
        auto tmp = gradient;
        gradient = 0;
        for (size_t j=0; j<gradient.N(); j++)
          for (size_t k=0; k<gradient.M(); k++)
            // Use sparseRange because jacobianTransposed may be a sparse DiagonalMatrix
            for(auto&& [jacobianInverseTransposed_j_l, l] : sparseRange(jacobianInverseTransposed[j]))
              gradient[j][k] += jacobianInverseTransposed_j_l * tmp[l][k];
      }
    }

    /** \brief Wrapper around a callable that applies the inverse Piola transform
     *
     * The LocalInterpolation implementations in dune-localfunctions expect local-valued
     * functions, but the ones dune-functions expect global-valued ones.  Therefore,
     * we need to stuff the inverse Piola transform between dune-functions and
     * dune-localfunctions, and this is what this class does.
     */
    template<class Function, class LocalCoordinate, class Element>
    class LocalValuedFunction
    {
      const Function& f_;
      const Element& element_;

    public:

      LocalValuedFunction(const Function& f, const Element& element)
      : f_(f), element_(element)
      {}

      auto operator()(const LocalCoordinate& xi) const
      {
        auto&& f = Dune::Impl::makeFunctionWithCallOperator<LocalCoordinate>(f_);
        auto globalValue = f(xi);

        // Apply the inverse Piola transform
        auto jacobianTransposed = element_.geometry().jacobianTransposed(xi);

        auto localValue = globalValue;
        jacobianTransposed.mv(globalValue, localValue);

        return localValue;
      }
    };
  };

  /** \brief Implementation of a dune-localfunctions LocalBasis that applies a transformation
   *
   * \tparam Transformator The transformation (e.g., Piola) that is to be applied
   * \tparam LocalValuedLocalBasis The local-valued LocalBasis this is getting transformed
   * \tparam Element The element that the global-valued FE lives on
   */
  template<class Transformator, class LocalValuedLocalBasis, class Element>
  class GlobalValuedLocalBasis
  {
  public:
    using Traits = typename LocalValuedLocalBasis::Traits;

    /** \brief Bind the local basis to a particular grid element
     */
    void bind(const LocalValuedLocalBasis& localValuedLocalBasis, const Element& element)
    {
      localValuedLocalBasis_ = &localValuedLocalBasis;
      element_ = &element;
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

      Transformator::apply(out, x, element_->geometry());
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

      Transformator::applyJacobian(out, x, element_->geometry());
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int,2>& order,
                 const typename Traits::DomainType& x,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(x, out);
      } else if (totalOrder == 1) {
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        out.resize(size());

        // TODO: The following is wasteful:  We compute the full Jacobian and then return
        // only a part of it.  While we need the full Jacobian of the underlying local-valued LFE,
        // it should be possible to compute only a partial Piola transform for the requested
        // partial derivatives.
        std::vector<typename Traits::JacobianType> fullJacobian;
        localValuedLocalBasis_->evaluateJacobian(x,fullJacobian);

        Transformator::applyJacobian(fullJacobian, x, element_->geometry());

        for (std::size_t i=0; i<out.size(); i++)
          for (std::size_t j=0; j<out[i].size(); j++)
            out[i][j] = fullJacobian[i][j][direction];

      } else
        DUNE_THROW(NotImplemented, "Partial derivatives of order 2 or higher");
    }

    //! \brief Polynomial order of the shape functions
    auto order() const
    {
      return localValuedLocalBasis_->order();
    }

    const LocalValuedLocalBasis* localValuedLocalBasis_;
    const Element* element_;
  };

  /** \brief Implementation of a dune-localfunctions LocalInterpolation
   *    that accepts global-valued functions
   *
   * \tparam Transformator The transformation (e.g., Piola) that transforms from local to global values
   * \tparam LocalValuedLocalInterpolation The local-valued LocalInterpolation
   *    that is used for the actual interpolation
   * \tparam Element The element that the global-valued FE lives on
   */
  template<class Transformator, class LocalValuedLocalInterpolation, class Element>
  class GlobalValuedLocalInterpolation
  {
  public:
    /** \brief Bind the local interpolation object to a particular grid element
     */
    void bind(const LocalValuedLocalInterpolation& localValuedLocalInterpolation, const Element& element)
    {
      localValuedLocalInterpolation_ = &localValuedLocalInterpolation;
      element_ = &element;
    }

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
      typename Transformator::template LocalValuedFunction<F,LocalCoordinate,Element> localValuedFunction(f, *element_);
      localValuedLocalInterpolation_->interpolate(localValuedFunction, out);
    }

  private:
    const LocalValuedLocalInterpolation* localValuedLocalInterpolation_;
    const Element* element_;
  };


  /** \brief LocalFiniteElement implementation that uses values defined wrt particular grid elements
   *
   * \tparam Transformator Class implementing range-space transformations (like the Piola transforms)
   * \tparam LocalValuedLFE LocalFiniteElement implementation whose values are to be transformed
   * \tparam Element Element where to transform the FE values to
   */
  template<class Transformator, class LocalValuedLFE, class Element>
  class GlobalValuedLocalFiniteElement
  {
    using LocalBasis = GlobalValuedLocalBasis<Transformator,
                                              typename LocalValuedLFE::Traits::LocalBasisType,
                                              Element>;
    using LocalInterpolation = GlobalValuedLocalInterpolation<Transformator,
                                                              typename LocalValuedLFE::Traits::LocalInterpolationType,
                                                              Element>;

  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<LocalBasis,
                                            typename LocalValuedLFE::Traits::LocalCoefficientsType,
                                            LocalInterpolation>;

    GlobalValuedLocalFiniteElement() {}

    void bind(const LocalValuedLFE& localValuedLFE, const Element& element)
    {
      globalValuedLocalBasis_.bind(localValuedLFE.localBasis(), element);
      globalValuedLocalInterpolation_.bind(localValuedLFE.localInterpolation(), element);
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
      return globalValuedLocalInterpolation_;
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

    typename Traits::LocalBasisType globalValuedLocalBasis_;
    typename Traits::LocalInterpolationType globalValuedLocalInterpolation_;
    const LocalValuedLFE* localValuedLFE_;
  };

}        // namespace Dune::Functions::Impl

#endif   // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GLOBALVALUEDLOCALFINITEELEMENT_HH
