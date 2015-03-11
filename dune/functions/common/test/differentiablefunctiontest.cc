// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <memory>
#include <functional>

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

#include <dune/functions/common/differentiablefunction.hh>

#include <dune/functions/analyticfunctions/polynomial.hh>
#include <dune/functions/analyticfunctions/trigonometricfunction.hh>

//#include <dune/functions/common/callable.hh>

//#include "derivativecheck.hh"


// Check if interface compiles and is implementable by a simple dummy
struct DifferentiableFunctionImplementableTest
{

  template<class F>
  static bool checkWithFunction(F&& f)
  {
    bool passed = true;

    {
      std::cout << "--------------" << std::endl;

//      passed = passed and DerivativeCheck<QuadraticPolynomial>::checkAllImplementedTrulyDerived(testFunction, 10);

      // Test whether I can evaluate the function somewhere
      std::cout << "Function value at x=5: " << f(5) << std::endl;



      std::cout << std::endl << "Check calling derivatives through function object" << std::endl;

      // Test whether I can evaluate the first derivative
      auto df = derivative(f);
      std::cout << "Derivative at x=5: " << df(5) << std::endl;

      // Test whether I can evaluate the second derivative through FunctionHandle
      auto ddf = derivative(df);
      std::cout << "Second derivative at x=5: " << ddf(5) << std::endl;

      // Test whether I can evaluate the third derivative through FunctionHandle
      auto dddf = derivative(ddf);
      std::cout << "Third derivative at x=5: " << dddf(5) << std::endl;



      std::cout << std::endl << "Check calling derivatives through DifferentiableFunction object" << std::endl;
      Dune::Functions::DifferentiableFunction<double(const double&)> fi = f;

      // Try to reassign wrapper
      fi = f;

      // Try assigning a default constructed wrapper
      Dune::Functions::DifferentiableFunction<double(const double&)> fii;
      fii = fi;

      // Try to copy wrapper
      auto fiii = fii;

      std::cout << "Function value at x=5: " << fiii(5) << std::endl;
      // Test whether I can evaluate the first derivative
      auto dfiii = derivative(fiii);
      std::cout << "Derivative at x=5: " << dfiii(5) << std::endl;

      // Test whether I can evaluate the second derivative through FunctionHandle
      auto ddfiii = derivative(dfiii);
      std::cout << "Second derivative at x=5: " << ddfiii(5) << std::endl;

      // Test whether I can evaluate the third derivative through FunctionHandle
      auto dddfiii = derivative(ddfiii);
      std::cout << "Third derivative at x=5: " << dddfiii(5) << std::endl;

      // Wrap as non-differentiable function
      Dune::Functions::DifferentiableFunction<double(const double&)> g = [=] (const double& x) {return f(x);};
      std::cout << "Function value at x=5: " << g(5) << std::endl;

      try {
        auto dg = derivative(g);
        passed = false;
      }
      catch (Dune::NotImplemented e)
      {
        std::cout << "Obtaining derivative from nondifferentiable function failed expectedly" << std::endl;
      }

    }

#if 0
    std::cout << std::endl << "Check calling persistent derivatives through shared_ptr" << std::endl;

    // Test whether I can evaluate the first derivative through shared_ptr
    persistentDerivative->evaluate(5, df);
    std::cout << "Derivative at x=5: " << df << std::endl;

    // Test whether I can evaluate the second derivative through shared_ptr
    persistentSecondDerivative->evaluate(5, ddf);
    std::cout << "Second derivative at x=5: " << ddf << std::endl;

    // Test whether I can evaluate the third derivative through shared_ptr
    persistentThirdDerivative->evaluate(5, dddf);
    std::cout << "Third derivative at x=5: " << dddf << std::endl;
#endif

    return passed;
  }

  static bool check()
  {
    bool passed = true;

    passed = passed and checkWithFunction(Dune::Functions::Polynomial<double>({1, 2, 3}));
    passed = passed and checkWithFunction(Dune::Functions::TrigonometricFunction<double, 1, 0>());

    return passed;
  }


};



#if 0
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
#endif





int main ( int argc, char **argv )
try
{
  bool passed = true;

  passed = passed and DifferentiableFunctionImplementableTest::check();

#if 0
  passed = passed and DerivativeRangeTerminationTest<double, double, 5>::check();

  passed = passed and DerivativeRangeTerminationTest<Dune::FieldVector<double, 3> , Dune::FieldVector<double, 1>, 5 >::check();

  passed = passed and DerivativeRangeTerminationTest<Dune::FieldVector<double, 1> , Dune::FieldVector<double, 1>, 5 >::check();
#endif


  if (passed)
    std::cout << "All tests passed" << std::endl;


  return passed ? 0: 1;
}
catch( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
