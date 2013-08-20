// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <memory>

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

#include <dune/functions/common/differentiablefunction.hh>


class QuadraticPolynomial
: public Dune::Functions::DifferentiableFunction<double,double>
{
public:

  /** \brief Constructor
   * \param a,b,c Coefficients with respect to the monomial basis
   */
  QuadraticPolynomial(double a, double b, double c)
  : a_(a), b_(b), c_(c)
  {}

  /** \brief Evaluate the function at a given point */
  void evaluate(const double& x, double& f) const
  {
    f = a_*x*x + b_*x + c_;
  }

  /** \brief Get the function implementing the first derivative */
  QuadraticPolynomial* derivative() const
  {
    if (not derivative_)
      derivative_ = Dune::make_shared<QuadraticPolynomial>(0, 2*a_, b_);
    return derivative_.get();
  }

private:
  // coefficients
  double a_, b_, c_;

  mutable Dune::shared_ptr<QuadraticPolynomial> derivative_;

};



int main ( int argc, char **argv )
try
{
  QuadraticPolynomial testFunction(1,1,1);

  // Test whether I can evaluate the function somewhere
  double f;
  testFunction.evaluate(5, f);
  std::cout << "Function value at x=5: " << f << std::endl;

  // Test whether I can evaluate the first derivative
  auto derivative = Dune::Functions::derivative(testFunction);
  double df;
  derivative.evaluate(5, df);
  std::cout << "Derivative at x=5: " << df << std::endl;

  // Test whether I can evaluate the first derivative
  auto secondDerivative = Dune::Functions::derivative(derivative);
  double ddf;
  secondDerivative.evaluate(5, ddf);
  std::cout << "Second derivative at x=5: " << ddf << std::endl;

  return 0;
}
catch( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
}
