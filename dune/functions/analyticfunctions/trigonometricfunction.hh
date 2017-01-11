// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_ANALYTICFUNCTIONS_TRIGONOMETRICFUNCTION_HH
#define DUNE_FUNCTIONS_ANALYTICFUNCTIONS_TRIGONOMETRICFUNCTION_HH


namespace Dune {
namespace Functions {



/**
 * \brief A linear combination of trigonomic functions
 *
 * \ingroup FunctionImplementations
 *
 * \tparam K Scalar type. The polynomial will map K to K
 * \tparam sinFactor Factor in front of sin
 * \tparam cosFactor Factor in front of cos
 *
 * This class exists mainly to demonstrate how to implement
 * the \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> concept.
 */
template<class K, int sinFactor, int cosFactor>
class TrigonometricFunction
{
public:
  //! Evaluate function
  K operator () (const K& x) const
  {
    return sinFactor * std::sin(x) + cosFactor * std::cos(x);
  }
};


//! Obtain derivative of TrigonometricFunction function \ingroup FunctionImplementations
template<class K, int sinFactor, int cosFactor>
TrigonometricFunction<K, -cosFactor, sinFactor> derivative(const TrigonometricFunction<K, sinFactor, cosFactor>& f)
{
  return TrigonometricFunction<K, -cosFactor, sinFactor>();
}



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_ANALYTICFUNCTIONS_TRIGONOMETRICFUNCTION_HH
