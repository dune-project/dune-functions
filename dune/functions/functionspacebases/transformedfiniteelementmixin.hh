// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDFINITEELEMENTMIXIN_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDFINITEELEMENTMIXIN_HH

namespace Dune::Functions::Impl{
/** \brief Implementation of a dune-localfunctions LocalBasis that applies a
 * linear transformation
 *
 * \tparam FE The Finite Element
 */
template<class FE, class ReferenceLocalBasisTraits>
class TransformedLocalBasis
{

  public:
    using Traits = ReferenceLocalBasisTraits;

    TransformedLocalBasis(FE const& fe_impl) : fe_impl_(&fe_impl) {}

  public:
    /** \brief Number of shape functions
        This does generally not equal the size of the underlying local finite
       element.
     */
    auto size() const
    {
      return fe_impl_->size();
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType &x,
                          std::vector<typename Traits::RangeType> &out) const
    {

      rangeBuffer_.resize(fe_impl_->size());
      fe_impl_->referenceLocalBasis().evaluateFunction(x, rangeBuffer_);
      out.resize(size());
      fe_impl_->transform(rangeBuffer_, out);
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference element where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType &x,
                          std::vector<typename Traits::JacobianType> &out) const
    {

      jacobianBuffer_.resize(fe_impl_->size());
      fe_impl_->referenceLocalBasis().evaluateJacobian(x, jacobianBuffer_);
      out.resize(size());
      fe_impl_->transform(jacobianBuffer_, out);
    }

    /** \brief Evaluate Hessian of all shape functions
     *   \note Sfinae protected: this method is only available, if the wrapped
     * LocalBasis 1. exports a <HessianType> and 2. provides a evaluateHessian
     * method with a corresponding signature \param x Point in the reference
     * element where to evaluation the Hessians \param[out] out The Hessians of
     * all shape functions at the point x
     */
    template<class TT,
             std::enable_if_t<std::is_same_v<TT, typename Traits::HessianType>, int> = 0>
    void evaluateHessian(const typename Traits::DomainType &x,
                         std::vector<TT> &out) const
    {
      hessianBuffer_.resize(fe_impl_->size());
      fe_impl_->referenceLocalBasis().evaluateHessian(x, hessianBuffer_);
      out.resize(size());
      fe_impl_->transform(hessianBuffer_, out);
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index
     * notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(std::array<unsigned int, Traits::dimDomain> const &order,
                 const typename Traits::DomainType &x,
                 std::vector<typename Traits::RangeType> &out) const
    {

      rangeBuffer_.resize(fe_impl_->size());
      fe_impl_->referenceLocalBasis().partial(order, x, rangeBuffer_);
      out.resize(size());
      fe_impl_->transform(rangeBuffer_, out);
    }

    //! \brief Polynomial order of the shape functions
    auto order() const { return fe_impl_->referenceLocalBasis().order(); }

    //Transformator const &transformator() const { return *transformator_; }

    // auto const &cacheable() const { return fe_impl_->referenceBasis(); }

  private:
    FE const* fe_impl_;
    mutable std::vector<typename Traits::RangeType> rangeBuffer_;
    mutable std::vector<typename Traits::JacobianType> jacobianBuffer_;
    mutable std::vector<typename Traits::HessianType> hessianBuffer_;
};


template<class FE, class ReferenceLocalBasisTraits>
class TransformedFiniteElementMixin
{
  public:
    TransformedFiniteElementMixin()
    : tlb_(this->as_impl()){}

    TransformedFiniteElementMixin(TransformedFiniteElementMixin const& other)
    : TransformedFiniteElementMixin(){}

    FE  const& as_impl() const { return *(static_cast<FE const*>(this));}

    auto const& localBasis() const{ return tlb_;}

  protected:
    TransformedLocalBasis<FE, ReferenceLocalBasisTraits> tlb_;
};


// Example Usage
template<class ReferenceLocalBasisTraits>
class TransformedFEExample : TransformedFiniteElementMixin<TransformedFEExample<ReferenceLocalBasisTraits>, ReferenceLocalBasisTraits>
{
    friend class TransformedLocalBasis<TransformedFEExample, ReferenceLocalBasisTraits>;
    using ReferenceRangeType = double;
    using RangeType = double;
    protected:
        auto const& referenceLocalBasis(){/* return the basis that is to be transformed*/}
        void transform(std::vector<ReferenceRangeType> const& in, std::vector<RangeType>& out){
            // apply transformation
        }
    public:
        auto localCoefficients(){ /* return local Coefficients of the transformed FE*/}
        auto localInterpolation(){ /* return local Interpolation fitting to the transformed FE(, i.e. the push forwards of the global DoFs onto the reference Element*/}
        auto size(){ /* return the size of the finite element after transformation)*/}

        template<class Element>
        void bind(Element const& e){
            // calculate transformation T
            // get local state if needed
        }


};
} // namespace Dune::Functions::Impl
#endif