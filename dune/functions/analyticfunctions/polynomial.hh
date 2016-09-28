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
 * \tparam K Scalar type. The polynomial will map K to K
 *
 * This class exists mainly to demonstrate how to implement
 * the \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> concept.
 */
template<class K>
class Polynomial
{
public:

  //! Default constructor
  Polynomial() = default;

  //! Copy constructor
  Polynomial(const Polynomial& other) = default;

  //! Move constructor
  Polynomial(Polynomial&& other) = default;

  /**
   * \brief Create from list of coefficients
   *
   * Coefficients are ordered in accordance with
   * the corresponding monomial order
   */
  Polynomial(std::initializer_list<double> coefficients) :
    coefficients_(coefficients)
  {}

  /**
   * \brief Create from list of coefficients
   *
   * Coefficients are ordered in accordance with
   * the corresponding monomial order. This will
   * move the coefficients from givven vector.
   */
  Polynomial(std::vector<K>&& coefficients) :
      coefficients_(std::move(coefficients))
  {}

  /**
   * \brief Create from list of coefficients
   *
   * Coefficients are ordered in accordance with
   * the corresponding monomial order. This will
   * copy coefficients from given vector.
   */
  Polynomial(const std::vector<K>& coefficients) :
    coefficients_(coefficients)
  {}

  //! Evaluate polynomial
  K operator() (const K& x) const
  {
    auto y = K(0);
    for (size_t i=0; i<coefficients_.size(); ++i)
      y += coefficients_[i] * std::pow(x, i);
    return y;
  }

  /**
   * \brief Obtain derivative of Polynomial function
   *
   * \ingroup FunctionImplementations
   *
   * The derivative contains its own coefficient
   * list and is not updated if the original function
   * is changed.
   */
  friend Polynomial derivative(const Polynomial& p)
  {
    std::vector<K> dpCoefficients(p.coefficients().size()-1);
    for (size_t i=1; i<p.coefficients_.size(); ++i)
      dpCoefficients[i-1] = p.coefficients()[i]*i;
    return Polynomial(std::move(dpCoefficients));
  }

  //! Obtain reference to coefficient vector
  const std::vector<K>& coefficients() const
  {
    return coefficients_;
  }

private:
  std::vector<K> coefficients_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_ANALYTICFUNCTIONS_POLYNOMIAL_HH
