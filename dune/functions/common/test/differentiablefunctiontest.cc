// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <memory>
//#include <functional>

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/callable.hh>

#include "derivativecheck.hh"

// Check if interface compiles and is implementable by a simple dummy
struct DifferentiableFunctionImplementableTest
{

  class QuadraticPolynomial
  : public Dune::Functions::DifferentiableFunction<double,double>
  {
  public:

    // Important: Explicitly export exact derivative type
    typedef QuadraticPolynomial Derivative;

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

  static bool check()
  {
    bool passed = true;

    QuadraticPolynomial testFunction(1,1,1);

    // Test whether I can evaluate the function somewhere
    double f;
    testFunction.evaluate(5, f);
    std::cout << "Function value at x=5: " << f << std::endl;



    std::cout << std::endl << "Check calling derivatives through FunctionHandle" << std::endl;

    // Test whether I can evaluate the first derivative
    auto derivative = Dune::Functions::derivative(testFunction);
    double df;
    derivative.evaluate(5, df);
    std::cout << "Derivative at x=5: " << df << std::endl;

    // Test whether I can evaluate the second derivative through FunctionHandle
    auto secondDerivative = Dune::Functions::derivative(derivative);
    double ddf;
    secondDerivative.evaluate(5, ddf);
    std::cout << "Second derivative at x=5: " << ddf << std::endl;

    // Test whether I can evaluate the third derivative through FunctionHandle
    auto thirdDerivative = Dune::Functions::derivative(secondDerivative);
    double dddf;
    thirdDerivative.evaluate(5, dddf);
    std::cout << "Third derivative at x=5: " << dddf << std::endl;



    std::cout << std::endl << "Check calling persistent derivatives through shared_ptr" << std::endl;

    // Test whether I can evaluate the first derivative through shared_ptr
    auto persistentDerivative = derivative.shared_ptr();
    persistentDerivative->evaluate(5, df);
    std::cout << "Derivative at x=5: " << df << std::endl;

    // Test whether I can evaluate the second derivative through shared_ptr
    auto persistentSecondDerivative = secondDerivative.shared_ptr();
    persistentSecondDerivative->evaluate(5, ddf);
    std::cout << "Second derivative at x=5: " << ddf << std::endl;

    // Test whether I can evaluate the third derivative through shared_ptr
    auto persistentThirdDerivative = thirdDerivative.shared_ptr();
    persistentThirdDerivative->evaluate(5, dddf);
    std::cout << "Third derivative at x=5: " << dddf << std::endl;



    std::cout << std::endl << "Check calling function and derivatives through Callable wrapper" << std::endl;

    auto callableF = Dune::Functions::callable(testFunction);
    std::cout << "Function value at x=5: " << callableF(5) << std::endl;

    auto callableDF = Dune::Functions::callable(Dune::Functions::derivative(testFunction));
    std::cout << "Derivative at x=5: " << callableDF(5) << std::endl;

    auto callableDDF = Dune::Functions::callable(Dune::Functions::derivative(Dune::Functions::derivative(testFunction)));
    std::cout << "Second derivative at x=5: " << callableDDF(5) << std::endl;



//    std::cout << std::endl << "Check calling function and derivatives through std::function" << std::endl;

//    std::function<double(double)> stdF = Dune::Functions::callable(testFunction);
//    std::cout << "Function value at x=5: " << stdF(5) << std::endl;

//    std::function<double(double)> stdDF = Dune::Functions::callable(Dune::Functions::derivative(testFunction));
//    std::cout << "Derivative at x=5: " << stdDF(5) << std::endl;

    passed = passed and DerivativeCheck<QuadraticPolynomial>::checkAllImplementatedTrulyDerived(testFunction, 10);

    return passed;
  }
};



// Check if recursive DerivativeRange definition terminates
// after at most maxRecursionLevel=k derivatives, i.e. if
// the type of the k-th and (k+1)-nd derivative is the same.
template<class D, class R, int maxRecursionLevel>
struct DerivativeRangeTerminationTest
{

  template<class DR, int recursionLevel>
  struct TerminationTest
  {
    static int level()
    {
      std::cout << "Type of " << recursionLevel << "-th derivative is " << Dune::className<DR>() << std::endl;
      typedef typename Dune::Functions::DerivativeTraits<D, DR>::DerivativeRange DDR;
      int upperBound = TerminationTest<DDR, recursionLevel+1>::level();
      if (Dune::is_same<DR, DDR>::value)
        return recursionLevel;
      else
        return upperBound;
    }
  };

  template<class DR>
  struct TerminationTest<DR, maxRecursionLevel>
  {
    static int level()
    {
      std::cout << "Type of " << maxRecursionLevel << "-th derivative is " << Dune::className<DR>() << std::endl;
      typedef typename Dune::Functions::DerivativeTraits<D, DR>::DerivativeRange DDR;
      if (Dune::is_same<DR, DDR>::value)
        return maxRecursionLevel;
      else
        return maxRecursionLevel+1;
    }
  };

  static bool check()
  {
    std::cout << "Checking recursion for Domain=" << Dune::className<D>() << " and Range=" << Dune::className<R>() << std::endl;
    int terminationLevel = TerminationTest<R, 0>::level();
    if (terminationLevel <= maxRecursionLevel)
    {
      std::cout << "Recursion terminated after " << terminationLevel << "-th derivative" << std::endl;
      return true;
    }
    else
    {
      std::cout << "Recursion did not terminated after given maxRecursionLevel " << maxRecursionLevel;
      return false;
    }
  }

};





int main ( int argc, char **argv )
try
{
  bool passed = true;

  passed = passed and DifferentiableFunctionImplementableTest::check();

  passed = passed and DerivativeRangeTerminationTest<double, double, 5>::check();

  passed = passed and DerivativeRangeTerminationTest<Dune::FieldVector<double, 3> , Dune::FieldVector<double, 1>, 5 >::check();

  passed = passed and DerivativeRangeTerminationTest<Dune::FieldVector<double, 1> , Dune::FieldVector<double, 1>, 5 >::check();


  if (passed)
    std::cout << "All tests passed" << std::endl;


  return passed ? 0: 1;
}
catch( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
}
