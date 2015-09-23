// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_ANALYTICFUNCTIONS_TRIGONOMETRICFUNCTION_HH
#define DUNE_FUNCTIONS_ANALYTICFUNCTIONS_TRIGONOMETRICFUNCTION_HH


namespace Dune {
namespace Functions {



template<class K, int sinFactor, int cosFactor>
class TrigonometricFunction
{
public:
  K operator () (const K& x) const
  {
    return sinFactor * std::sin(x) + cosFactor * std::cos(x);
  }
};


template<class K, int sinFactor, int cosFactor>
TrigonometricFunction<K, -cosFactor, sinFactor> derivative(const TrigonometricFunction<K, sinFactor, cosFactor>& f)
{
  return TrigonometricFunction<K, -cosFactor, sinFactor>();
}



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_ANALYTICFUNCTIONS_TRIGONOMETRICFUNCTION_HH
