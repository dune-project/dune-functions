// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFEREENTIONABEFUNCTIONFROMCALLABLES_HH
#define DUNE_FUNCTIONS_COMMON_DIFFEREENTIONABEFUNCTIONFROMCALLABLES_HH


namespace Dune {
namespace Functions {



/**
 * \brief Wrap a list of callable objects as function modelling the DifferentiableFunction interface
 *
 * You can use this to implement a DifferentiableFunction including
 * a variable number of derivatives using callable objects.
 *
 * Note that using makeDifferentiableFunction will be less verbose than
 * creating this wrapper manually.
 */
template<class Signature, template<class> class DerivativeTraits, class... Callables>
class DifferentiableFunctionFromCallables;



template<class Range, class Domain, template<class> class DerivativeTraits, class F>
class DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F>
{
public:

  using Signature = Range(Domain);
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

  using Derivative = DifferentiableFunction<DerivativeSignature, DerivativeTraits>;

  template<class FF, disableCopyMove<DifferentiableFunctionFromCallables, FF> = 0>
  DifferentiableFunctionFromCallables(FF&& f) :
    f_(std::forward<FF>(f))
  {}

  Range operator() (const Domain& x) const
  {
    return f_(x);
  }

  friend Derivative derivative(const DifferentiableFunctionFromCallables& t)
  {
    DUNE_THROW(Dune::NotImplemented, "Derivative not implemented");
  }

private:
  F f_;
};



template<class Range, class Domain, template<class> class DerivativeTraits, class F, class DF, class... Derivatives>
class DifferentiableFunctionFromCallables<Range(Domain), DerivativeTraits, F, DF, Derivatives...>
{
public:

  using Signature = Range(Domain);
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

  using Derivative = DifferentiableFunctionFromCallables<DerivativeSignature, DerivativeTraits, DF, Derivatives...>;

  template<class FF, class DFF, class... DDFF>
  DifferentiableFunctionFromCallables(FF&& f, DFF&& df, DDFF&&... ddf) :
    f_(std::forward<FF>(f)),
    df_(std::forward<DFF>(df), std::forward<DDFF>(ddf)...)
  {}

  Range operator() (const Domain& x) const
  {
    return f_(x);
  }

  friend Derivative derivative(const DifferentiableFunctionFromCallables& t)
  {
    return t.df_;
  }

private:
  F f_;
  Derivative df_;
};


template<class Signature, template<class> class DerivativeTraits=DefaultDerivativeTraits>
struct SignatureTag;

/**
 * \brief Tag-class to encapsulate signature information
 *
 * \tparam Range range type
 * \tparam Domain domain type
 * \tparam DerivativeTraits traits template used to determine derivative traits
 */
template<class Range, class Domain, template<class> class DerivativeTraitsT>
struct SignatureTag<Range(Domain), DerivativeTraitsT>
{
  using Signature = Range(Domain);

  template<class T>
  using DerivativeTraits = DerivativeTraitsT<T>;
};


/**
 * \brief Create a DifferentiableFunction from callables
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
