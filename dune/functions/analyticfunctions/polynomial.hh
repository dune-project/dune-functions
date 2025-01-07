// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_ANALYTICFUNCTIONS_POLYNOMIAL_HH
#define DUNE_FUNCTIONS_ANALYTICFUNCTIONS_POLYNOMIAL_HH

#include <cmath>
#include <initializer_list>
#include <vector>


#include <dune/common/hybridutilities.hh>

namespace Dune {
namespace Functions {

namespace Impl {

  // Compute coefficients of derivative of polynomial.
  // Overload for std::vector
  template<class K, class Allocator>
  auto polynomialDerivativeCoefficients(const std::vector<K, Allocator>& coefficients) {
    if (coefficients.size()==0)
      return std::vector<K, Allocator>();
    std::vector<K, Allocator> dpCoefficients(coefficients.size()-1);
    for (size_t i=1; i<coefficients.size(); ++i)
      dpCoefficients[i-1] = coefficients[i]*K(i);
    return dpCoefficients;
  }

  // Compute coefficients of derivative of polynomial.
  // Overload for std::array
  template<class K, std::size_t n>
  auto polynomialDerivativeCoefficients(const std::array<K, n>& coefficients) {
    if constexpr (n==0)
      return coefficients;
    else
    {
      std::array<K, n-1> dpCoefficients;
      for (size_t i=1; i<coefficients.size(); ++i)
        dpCoefficients[i-1] = coefficients[i]*K(i);
      return dpCoefficients;
    }
  }

  // Compute coefficients of derivative of polynomial.
  // Helper function for the std::integer_sequence overload.
  // With C++20 this can be avoided, because lambda function
  // can partially specify template arguments which allows
  // to do the same inline.
  template<class I, I i0, I... i, class J, J j0, J... j>
  auto polynomialDerivativeCoefficientsHelper(std::integer_sequence<I, i0, i...>, std::integer_sequence<J, j0, j...>) {
    return std::integer_sequence<I, i*I(j)...>();
  }

  // Compute coefficients of derivative of polynomial.
  // Overload for std::integer_sequence
  template<class I, I... i>
  auto polynomialDerivativeCoefficients(std::integer_sequence<I, i...> coefficients) {
    if constexpr (sizeof...(i)==0)
      return coefficients;
    else
      return polynomialDerivativeCoefficientsHelper(coefficients, std::make_index_sequence<sizeof...(i)>());
  }

  // Compute coefficients of derivative of polynomial.
  // Overload for std::tuple
  template<class...T>
  auto polynomialDerivativeCoefficients(const std::tuple<T...>& coefficients) {
    if constexpr (sizeof...(T)==0)
      return coefficients;
    else
    {
      // Notice that std::multiplies<void> has issues with signed types.
      // E.g., `decltype(-2,2ul)` is `long unsigned int`.
      // Hence the same is deduced as return type in std::multiplies.
      // To avoid this, we explicitly pass the exponent `i+1` as signed type.
      // If the coefficient is signed, both types are now signed and
      // so is the deduced result type of std::multiplies.
      auto mult = Dune::Hybrid::hybridFunctor(std::multiplies());
      return Dune::unpackIntegerSequence([&](auto... i) {
        return std::tuple(mult(std::get<i+1>(coefficients),
          std::integral_constant<long signed int, (long signed int)(i+1)>()) ...);
      }, std::make_index_sequence<sizeof...(T)-1>());
    }
  }

} // namespace Impl in Dune::Functions::



/**
 * \brief A univariate polynomial implementation
 *
 * \ingroup FunctionImplementations
 *
 * \tparam K Scalar type. The polynomial will map K to K
 * \tparam C Coefficient container type (default std::vector<K>)
 *
 * This class will store a coefficient container of type C.
 * Supported containers are std::vector, std::array, std::integer_sequence, std::tuple.
 * When passing std::tuple, coefficients of type std::integral_constant are promoted
 * as std::integral_constant when computing derivatives.
 *
 * Class template argument deduction is supported for passing std::array, std::vector,
 * std::integer_sequence, std::initializer_list. When passing such containers
 * without specifying the template parameters, then the scalar type is deduced to
 * be the coefficient type. Notice that the deduced coefficient container type when
 * passing std::initializer_list is std::vector.
 *
 * If you want to use different types for scalar and coefficients, you can use the
 * makePolynomial() function to explicitly specify the scalar type while the
 * coefficient type is deduced.
 *
 * This class exists mainly to demonstrate how to implement
 * the \ref Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> concept.
 */
template<class K, class C=std::vector<K>>
class Polynomial
{
  template<class CC>
  struct IsIntegerSequence : public std::false_type {};

  template<class I, I... i>
  struct IsIntegerSequence<std::integer_sequence<I, i...>> : public std::true_type {};

public:

  //! The type of the stored coefficient container
  using Coefficients = C;

  //! Default constructor
  Polynomial() = default;

  /**
   * \brief Create from container of coefficients
   *
   * Coefficients are ordered in accordance with
   * the corresponding monomial order. The constructed
   * Polynomial object will store a copy of the passed
   * coefficient container.
   */
  Polynomial(Coefficients coefficients) :
      coefficients_(std::move(coefficients))
  {}

  //! Evaluate polynomial
  K operator() (const K& x) const
  {
    auto y = K(0);
    auto n = Dune::Hybrid::size(coefficients_);
    Dune::Hybrid::forEach(Dune::range(n), [&](auto i) {
      y += Dune::Hybrid::elementAt(coefficients_, i) * std::pow(x, int(i));
    });
    return y;
  }

  //! Comparison of coefficients
  bool operator==(const Polynomial& other) const
  {
    if constexpr (IsIntegerSequence<Coefficients>::value)
      return true;
    else
      return coefficients()==other.coefficients();
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
  friend auto derivative(const Polynomial& p)
  {
    auto derivativeCoefficients = Impl::polynomialDerivativeCoefficients(p.coefficients());
    using DerivativeCoefficients = decltype(derivativeCoefficients);
    return Polynomial<K, DerivativeCoefficients>(std::move(derivativeCoefficients));
  }

  //! Obtain reference to coefficient vector
  const Coefficients& coefficients() const
  {
    return coefficients_;
  }

private:
  Coefficients coefficients_;
};



template<class K>
Polynomial(std::vector<K>) -> Polynomial<K, std::vector<K>>;

template<class K, std::size_t n>
Polynomial(std::array<K,n>) -> Polynomial<K, std::array<K,n>>;

template<class K, K... ci>
Polynomial(std::integer_sequence<K, ci...>) -> Polynomial<K, std::integer_sequence<K,ci...>>;

template<class K>
Polynomial(std::initializer_list<K>) -> Polynomial<K, std::vector<K>>;



/**
 * \brief Create Polynomial
 *
 * \tparam K Scalar type. The polynomial will map K to K
 * \tparam C Coefficient container type
 *
 * This helper function allows to specify K, but lets C be deduced from
 * the argument. If the scalar type K is the same as the coefficient
 * type (i.e. the entry type of the coefficient container C), then
 * you can also rely on class template argument deduction and
 * call Polynomial(coefficients) directly.
 */
template<class K, class Coefficients>
auto makePolynomial(Coefficients coefficients)
{
  return Polynomial<K, Coefficients>(std::move(coefficients));
}

/**
 * \brief Create Polynomial
 *
 * \tparam K Scalar type. The polynomial will map K to K
 * \tparam C Coefficient type
 *
 * The initializer list will be stored as std::vector
 * in the created object of type Polynomial<K,std::vector<C>>.
 */
template<class K, class C>
auto makePolynomial(std::initializer_list<C> coefficients)
{
  return Polynomial<K>(std::move(coefficients));
}





}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_ANALYTICFUNCTIONS_POLYNOMIAL_HH
