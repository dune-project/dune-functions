// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFEREENTIONABEFUNCTIONFROMCALLABLES_HH
#define DUNE_FUNCTIONS_COMMON_DIFFEREENTIONABEFUNCTIONFROMCALLABLES_HH


#include <dune/common/typeutilities.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/functions/common/signature.hh>

#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/functionconcepts.hh>



namespace Dune {
namespace Functions {



template<class Signature, template<class> class DerivativeTraits, class... Callables>
class DifferentiableFunctionFromCallables;



/**
 * \brief Wrap a list of callable objects as derivative sequence modelling \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits>
 *
 * \ingroup FunctionImplementations
 *
 * \tparam Range Range type of function
 * \tparam Domain Domain type of function
 *
 * You can use this to implement a differentiable function including
 * a variable number of derivatives using callable objects.
 *
 * This models the \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> concept.
 *
 * Note that using makeDifferentiableFunction will be less verbose than
 * creating this wrapper manually.
 */
template<class Range, class Domain, template<class> class DerivativeTraits, class F>
class DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F>
{
public:

  //! Signature of function
  using Signature = Range(Domain);

  using RawSignature = typename SignatureTraits<Signature>::RawSignature;

  //! Signature of derivative
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

  //! Type of derivative
  using Derivative = DifferentiableFunction<DerivativeSignature, DerivativeTraits>;

  //! Constructor copying the given function
  template<class FF, disableCopyMove<DifferentiableFunctionFromCallables, FF> = 0>
  DifferentiableFunctionFromCallables(FF&& f) :
    f_(std::forward<FF>(f))
  {}

  //! Evaluate function
  Range operator() (const Domain& x) const
  {
    return f_(x);
  }

  /**
   * \brief Get derivative of DifferentiableFunctionFromCallables
   *
   * \ingroup FunctionImplementations
   */
  friend Derivative derivative(const DifferentiableFunctionFromCallables& t)
  {
    DUNE_THROW(Dune::NotImplemented, "Derivative not implemented");
  }

private:
  F f_;
};



/**
 * \brief Wrap a list of callable objects as derivative sequence modelling \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits>
 *
 * \ingroup FunctionImplementations
 *
 * \tparam Range Range type of function
 * \tparam Domain Domain type of function
 *
 * You can use this to implement a differentiable function including
 * a variable number of derivatives using callable objects.
 *
 * This models the \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> concept.
 *
 * Note that using makeDifferentiableFunction will be less verbose than
 * creating this wrapper manually.
 */
template<class Range, class Domain, template<class> class DerivativeTraits, class F, class DF, class... Derivatives>
class DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F, DF, Derivatives...>
{
public:

  using Signature = Range(Domain);
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

  using Derivative = DifferentiableFunctionFromCallables<DerivativeSignature, DerivativeTraits, DF, Derivatives...>;

  /**
   * \brief Constructor copying the given functions
   *
   * The arguments are used as implementation of the functions itself
   * and its derivatives with increasing order
   */
  template<class FF, class DFF, class... DDFF>
  DifferentiableFunctionFromCallables(FF&& f, DFF&& df, DDFF&&... ddf) :
    f_(std::forward<FF>(f)),
    df_(std::forward<DFF>(df), std::forward<DDFF>(ddf)...)
  {}

  //! Evaluate function
  Range operator() (const Domain& x) const
  {
    return f_(x);
  }

  /**
   * \brief Get derivative of DifferentiableFunctionFromCallables
   *
   * \ingroup FunctionImplementations
   */
  friend Derivative derivative(const DifferentiableFunctionFromCallables& t)
  {
    return t.df_;
  }

private:
  F f_;
  Derivative df_;
};


/**
 * \brief Create a DifferentiableFunction from callables
 *
 * \ingroup FunctionImplementations
 *
 * This will return a wrapper modelling the DifferentiableFunction interface
 * where the evaluation of the function and its derivatives are implemented
 * by the given callable objects.
 *
 * \param signatureTag A dummy parameter to pass the signature and derivative traits
 * \param f Callable objects implementing the evaluation of the function and its derivatives
 *
 * \returns Object modelling DifferentiableFunction interface
 */
template<class Signature, template<class> class DerivativeTraits, class... F>
DifferentiableFunctionFromCallables<Signature, DerivativeTraits, F...>
  makeDifferentiableFunctionFromCallables(const SignatureTag<Signature, DerivativeTraits>& signatureTag, F&&... f)
{
  return DifferentiableFunctionFromCallables<Signature, DerivativeTraits, F...>(f...);
}



} // namespace Functions
} // namespace Dune

#endif //DUNE_FUNCTIONS_COMMON_DIFFEREENTIONABEFUNCTIONFROMCALLABLES_HH
