// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BOUND_TEST_LOCALFE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BOUND_TEST_LOCALFE_HH

/** \file \brief Unit tests for GlobalValuedLocalFiniteElement objects that utilize
 * global(!)derivative degrees of freedom
 *
 *  \note This file is an adapted copy from
 * dune-localfunctions/dune/localfunctions/test/test-localfe.hh The tests in this file do not
 * replace the tests in the file mentioned above, but apply them on the local finite element bound
 * to a grid element. For the same reason this file is part of dune-functions rather than dune-
 * localfunctions.
 *
 *  \note This header is not part of the official Dune API and might be subject to
 * change.  You can use this header to test external finite element implementations, but be warned
 * that your tests might break with future Dune versions.
 */

#include <array>
#include <iomanip>
#include <iostream>
#include <typeinfo>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/transpose.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/functions/common/defaultderivativetraits.hh>

// Shapefunctions need to have a derivative
template <class LFE, class Element>
class ShapeFunctionDerivativeAsCallable;

// This class wraps one shape function of a local finite element as a callable with derivative() method
// that can be fed to the LocalInterpolation::interpolate method, even if this evaluates derivatives.
// implements Dune::Functions::DifferentiableFunction
template <class FE, class Element>
class ShapeFunctionAsCallableWithDerivative
{
public:
  using DomainType = typename FE::Traits::LocalBasisType::Traits::DomainType ;
  using RangeType = typename FE::Traits::LocalBasisType::Traits::RangeType ;
  using RangeFieldType = typename FE::Traits::LocalBasisType::Traits::RangeFieldType ;
  static constexpr int dim = FE::Traits::LocalBasisType::Traits::dimDomain;
  static int constexpr dimRange = FE::Traits::LocalBasisType::Traits::dimRange;

  ShapeFunctionAsCallableWithDerivative(const FE &fe, int shapeFunction, Element const &element)
    : fe_(fe)
    , shapeFunction_(shapeFunction)
    , e_(element)
  {}

  RangeType operator()(DomainType x) const
  {
    std::vector<RangeType> yy;
    fe_.localBasis().evaluateFunction(x, yy);
    return yy[shapeFunction_];
  }

  friend ShapeFunctionDerivativeAsCallable<FE, Element> derivative(ShapeFunctionAsCallableWithDerivative const &t)
  {
    return ShapeFunctionDerivativeAsCallable(t);
  }

public:
  const FE &fe_;
  int shapeFunction_;
  Element const &e_;
};

// Test whether a Finite Element has an `evaluateHessian(...)` method and exports `HessianType`.
template<class FE, class = void>
struct hasEvaluateHessian{
  static constexpr bool value = false;
};

template<class FE>
struct hasEvaluateHessian<FE, std::void_t<decltype(std::declval<FE>().evaluateHessian(FE::Traits::LocalBasisType::Traits::DomainType, FE::Traits::LocalBasisType::Traits::HessianType))>>
{
  static constexpr bool value = true;
};

// The derivative of one shape function
// implements Dune::Functions::DifferentiableFunction
template<class FE, class Element>
class ShapeFunctionDerivativeAsCallable : public ShapeFunctionAsCallableWithDerivative<FE, Element>
{

  using Base = ShapeFunctionAsCallableWithDerivative<FE, Element>;
  typedef typename FE::Traits::LocalBasisType::Traits::JacobianType JacobianType;

public:
  ShapeFunctionDerivativeAsCallable(Base const& base)
  :Base(base)
  {
  }

  JacobianType operator()(typename Base::DomainType const &x) const
  {
    std::vector<JacobianType> yy;
    Base::fe_.localBasis().evaluateJacobian(x, yy);
    return yy[Base::shapeFunction_] * Base::e_.geometry().jacobianInverse(x);
  }

  auto friend derivative(ShapeFunctionDerivativeAsCallable const &t)
  {
    if constexpr (hasEvaluateHessian<FE>::value){
      return [t](typename Base::DomainType const& x){

        std::array<typename FE::Traits::LocalBasisType::Traits::HessianType, Base::dimRange> hessian;
        std::vector<typename FE::Traits::LocalBasisType::Traits::HessianType> referenceHessians;
        t.fe_.localBasis().evaluateHessian(x, referenceHessians);
        const auto geometryJacobianInverse = t.e_.geometry().jacobianInverse(x);
        for (std::size_t k = 0; k < Base::dimRange; ++k)
          hessian[k] = transpose(geometryJacobianInverse) * referenceHessians[k] * geometryJacobianInverse;
        return hessian;
      };
    }
    else {
      return [t](typename Base::DomainType const &x)
      {
        std::array<Dune::FieldMatrix<typename Base::RangeFieldType, Base::dim, Base::dim>, Base::dimRange> referenceHessian, hessian;
        std::vector<typename Base::RangeType> partials;
        std::array<unsigned int, Base::dim> dir;
        for (std::size_t i = 0; i < Base::dim; ++i)
          for (std::size_t j = 0; j < Base::dim; ++j)
          {
            dir = {};
            dir[i]++;
            dir[j]++;
            t.fe_.localBasis().partial(dir, x, partials);
            for (std::size_t k = 0; k < Base::dimRange; ++k)
            {
              hessian[k][i][j] = 0.0;
              referenceHessian[k][i][j] = partials[t.shapeFunction_][k];
            }
          }

        const auto geometryJacobianInverse = t.e_.geometry().jacobianInverse(x);
        for (std::size_t k = 0; k < Base::dimRange; ++k)
          hessian[k] = transpose(geometryJacobianInverse) * referenceHessian[k] * geometryJacobianInverse;
        return hessian;
      };
    }
  }

};
// Check whether the degrees of freedom computed by LocalInterpolation
// are dual to the shape functions.  See Ciarlet, "The Finite Element Method
// for Elliptic Problems", 1978, for details.
template <class FE, class Element>
bool testLocalInterpolation(const FE &fe, Element const &element)
{
  bool ret = true;
  std::vector<typename ShapeFunctionAsCallableWithDerivative<FE, Element>::RangeFieldType> coeff;
  for (size_t i = 0; i < fe.size(); ++i)
  {
    //////////////////////////////////////////////////////////////////////////////
    // Feed the shape functions to the 'interpolate' method in form in form of a callable.
    //////////////////////////////////////////////////////////////////////////////

    // The i-th shape function as a function that 'interpolate' can deal with
    ShapeFunctionAsCallableWithDerivative<FE, Element> sfAsCallable(fe, i, element);

    // Compute degrees of freedom for that shape function
    // We expect the result to be the i-th unit vector
    fe.localInterpolation().interpolate(sfAsCallable, coeff);

    // Check size of weight vector
    if (coeff.size() != fe.localBasis().size())
    {
      std::cout << "Bug in LocalInterpolation for finite element type " << Dune::className(fe)
                << std::endl;
      std::cout << "    Interpolation produces " << coeff.size() << " degrees of freedom"
                << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size() << std::endl;
      std::cout << std::endl;
      return false;
    }

    // Check if interpolation weights are equal to coefficients
    for (std::size_t j = 0; j < coeff.size(); ++j)
    {
      using std::abs;
      if (abs(coeff[j] - (i == j)) > TOL)
      {
        std::cout << std::setprecision(16);
        std::cout << "Bug in LocalInterpolation for finite element type " << Dune::className(fe)
                  << std::endl;
        std::cout << "    Degree of freedom " << j << " applied to shape function " << i
                  << " yields value " << coeff[j] << ", not the expected value " << (i == j)
                  << std::endl;
        std::cout << std::endl;
        ret = false;
      }
    }
  }
  return ret;
}

// Function representing a constant function
// implements Dune::Functions::DifferentiableFunction
template <class Domain, class Range>
class Constant
{
public:
  Constant(double val) : value(val) {}

  Range operator()(Domain const &x) const { return value; }

  friend Constant<Domain, typename Dune::Functions::DefaultDerivativeTraits<Range(Domain)>::Range>
  derivative(Constant const &t)
  {
    return {0.};
  }

private:
  Range value;
};

// Check whether the space spanned by the shape functions
// contains the constant functions
template <class FE>
bool testCanRepresentDifferentiableConstants(const FE &fe, unsigned order = 5)
{
  typedef typename FE::Traits::LocalBasisType LB;
  using RangeType = typename LB::Traits::RangeType;
  bool success = true;

  // Construct the constant '1' function
  Constant<typename LB::Traits::DomainType, RangeType> constantOne(1.);
  // Project the constant function onto the FE space
  std::vector<double> coefficients;
  fe.localInterpolation().interpolate(constantOne, coefficients);

  // A set of test points
  const auto &quad = Dune::QuadratureRules<double, LB::Traits::dimDomain>::rule(fe.type(), order);

  // Loop over all quadrature points
  for (size_t i = 0; i < quad.size(); i++)
  {

    // Get a test point
    const auto &testPoint = quad[i].position();

    // Compute value of the representation of constantOne at the test point
    std::vector<RangeType> values;
    fe.localBasis().evaluateFunction(testPoint, values);

    RangeType sum(0);
    for (size_t j = 0; j < values.size(); j++)
      sum += coefficients[j] * values[j];

    if ((RangeType(1.0)-sum).two_norm() > TOL)
    {
      std::cout << "Finite element type " << Dune::className(fe)
                << " cannot represent constant functions!" << std::endl;
      std::cout << "    At position: " << testPoint << "," << std::endl;
      std::cout << "    discrete approximation of the '1' function has value " << sum << std::endl;
      std::cout << std::endl;
      success = false;
    }

  } // Loop over all quadrature points

  return success;
}

/** \brief Call tests for given finite element on a grid element
 *  \relates testFE
 * \param derivativePointSkip This is a small predicate class that allows to skip certain
 *   points when testing the derivative implementations.  It exists because some
 *   finite elements are not everywhere differentiable, but we still want to run
 *   the tests for derivatives.  Rather than constructing special sets of test
 *   points that avoid the problematic parts of the domain, we simply skip
 *   all test points that happen to be somewhere where the shape functions are
 *   not differentiable.
 */
template <unsigned int diffOrder = 0, class FE, class Element>
bool testBoundFE(const FE &fe, Element const &element, char disabledTests = DisableNone,
    const std::function<bool(const typename FE::Traits::LocalBasisType::Traits::DomainType &)>
        derivativePointSkip
    = nullptr)
{
  // Order of the quadrature rule used to generate test points
  unsigned int quadOrder = 2;

  bool success = true;

  if (FE::Traits::LocalBasisType::Traits::dimDomain != fe.type().dim())
  {
    std::cout << "Bug in type() for finite element type " << Dune::className(fe) << std::endl;
    std::cout << "    Coordinate dimension is " << FE::Traits::LocalBasisType::Traits::dimDomain
              << std::endl;
    std::cout << "    but GeometryType is " << fe.type() << " with dimension " << fe.type().dim()
              << std::endl;
    success = false;
  }

  if (fe.size() != fe.localBasis().size())
  {
    std::cout << "Bug in finite element type " << Dune::className(fe) << std::endl;
    std::cout << "    Size reported by LocalFiniteElement is " << fe.size() << std::endl;
    std::cout << "    but size reported by LocalBasis is " << fe.localBasis().size() << std::endl;
    success = false;
  }

  // Make sure evaluateFunction returns the correct number of values
  std::vector<typename FE::Traits::LocalBasisType::Traits::RangeType> values;
  fe.localBasis().evaluateFunction(
      Dune::ReferenceElements<double, FE::Traits::LocalBasisType::Traits::dimDomain>::general(
          fe.type())
          .position(0, 0),
      values);

  if (values.size() != fe.size())
  {
    std::cout << "Bug in finite element type " << Dune::className(fe) << std::endl;
    std::cout << "    LocalFiniteElement.size() returns " << fe.size() << "," << std::endl;
    std::cout << "    but LocalBasis::evaluateFunction returns " << values.size() << " values!"
              << std::endl;
    success = false;
  }

  if (fe.size() != fe.localCoefficients().size())
  {
    std::cout << "Bug in finite element type " << Dune::className(fe) << std::endl;
    std::cout << "    Size reported by LocalFiniteElement is " << fe.size() << std::endl;
    std::cout << "    but size reported by LocalCoefficients is " << fe.localCoefficients().size()
              << std::endl;
    success = false;
  }

  const auto &lc = fe.localCoefficients();
  for (size_t i = 0; i < lc.size(); i++)
  {
    const auto &lk = lc.localKey(i);
    if (lk.codim() > fe.type().dim())
    {
      std::cout << "Bug in finite element type " << Dune::className(fe) << std::endl;
      std::cout << "    Codimension reported by localKey(" << i << ") is " << lk.codim()
                << std::endl;
      std::cout << "    but geometry is " << fe.type() << " with dimension " << fe.type().dim()
                << std::endl;
      success = false;
    }
  }

  if (not(disabledTests & DisableLocalInterpolation))
  {
    success = testLocalInterpolation<FE>(fe, element) and success;
  }

  if (not (disabledTests & DisableRepresentConstants))
  {
    success = testCanRepresentDifferentiableConstants<FE>(fe) and success;
  }

  if (not(disabledTests & DisableJacobian))
  {
    success = testJacobian<FE>(fe, quadOrder, derivativePointSkip) and success;
  }
  else
  {
    // make sure diffOrder is 0
    success = (diffOrder == 0) and success;
  }

  if (not(disabledTests & DisableEvaluate))
  {
    success = TestPartial::test<FE>(fe, TOL, jacobianTOL, diffOrder, quadOrder, derivativePointSkip) and success;
  }
  return success;
}

#endif
