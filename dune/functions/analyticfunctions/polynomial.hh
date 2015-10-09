// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_ANALYTICFUNCTIONS_POLYNOMIAL_HH
#define DUNE_FUNCTIONS_ANALYTICFUNCTIONS_POLYNOMIAL_HH


namespace Dune {
namespace Functions {


/**
 * \brief A scalar polynomial implementation
 *
 * \ingroup FunctionImplementations
 *
 * This class exists mainly to demonstrate how to implement
 * the \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> concept.
 */
template<class K>
class Polynomial
{
public:

  Polynomial()
  {}

  Polynomial(const Polynomial& other) :
    coefficients_(other.coefficients_)
  {}

  Polynomial(Polynomial&& other) :
    coefficients_(std::move(other.coefficients_))
  {}

  Polynomial(std::initializer_list<double> coefficients) :
    coefficients_(coefficients)
  {}

  Polynomial(std::vector<K>&& coefficients) :
      coefficients_(std::move(coefficients))
  {}

  Polynomial(const std::vector<K>& coefficients) :
    coefficients_(coefficients)
  {}

  K operator() (const K& x) const
  {
    auto y = K(0);
    for (size_t i=0; i<coefficients_.size(); ++i)
      y += coefficients_[i] * std::pow(x, i);
    return y;
  }

  friend Polynomial derivative(const Polynomial& p)
  {
    auto derivative = Polynomial();
    derivative.coefficients_.resize(p.coefficients_.size()-1);
    for (size_t i=1; i<p.coefficients_.size(); ++i)
      derivative.coefficients_[i-1] = p.coefficients_[i]*i;
    return derivative;
  }

  const std::vector<K>& coefficients() const
  {
    return coefficients_;
  }

private:
  std::vector<K> coefficients_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_ANALYTICFUNCTIONS_POLYNOMIAL_HH
