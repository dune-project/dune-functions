// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <memory>

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/functionfromcallable.hh>

#include "derivativecheck.hh"

struct Sinus
{
  double operator() (double x) const
  {
    return std::sin(x);
  }
};

double cosinus(double x)
{
  return std::cos(x);
}




template<class T, int n>
T energyNorm(const Dune::FieldMatrix<T, n, n>& M, const Dune::FieldVector<T, n>& x)
{
  Dune::FieldVector<T, n> y;
  M.mv(x, y);
  return 0.5 * x.dot(y);
}


// Check if interface compiles and is implementable by a simple dummy
struct FunctionFromCallableTest
{

  static bool check()
  {
    bool passed = true;
    {
      Sinus sinus;

      Dune::Functions::FunctionFromCallable<double, double> f(
          sinus,
          cosinus,
          [] (double x) -> double {return -std::sin(x);});

      passed = passed and DerivativeCheck<Dune::Functions::FunctionFromCallable<double, double> >::checkAllImplementedTrulyDerived(f, 10);

      double x=3.14159;
      double y;

      f.evaluate(x,y);
      std::cout << "Value of f at x=" << x <<" is " << y << std::endl;

      auto df = Dune::Functions::derivative(f);
      df->evaluate(x,y);
      std::cout << "Value of f' at x=" << x <<" is " << y << std::endl;

      auto ddf = Dune::Functions::derivative(df);
      ddf->evaluate(x,y);
      std::cout << "Value of f'' at x=" << x <<" is " << y << std::endl;
    }

    std::cout << std::endl;
    {
      typedef Dune::FieldVector<double, 3> Domain;
      typedef double Range;
      typedef Dune::FieldVector<double, 3> DRange;
      typedef Dune::FieldMatrix<double, 3, 3> DDRange;

      Dune::FieldMatrix<double,3,3> M(1.0);
      for(int i=0; i<M.N(); ++i)
        M[i][i] = 2.0;

      Dune::Functions::FunctionFromCallable<Domain, Range> f(
          std::bind(energyNorm<double, 3>, M, std::placeholders::_1),
          [=] (const Domain& x) -> DRange {
            DRange y;
            M.mv(x,y);
            return y;
          },
          [=] (const Domain& x) -> DDRange {
            return M;
          });

      passed = passed and DerivativeCheck<Dune::Functions::FunctionFromCallable<Domain, Range> >::checkAllImplementedTrulyDerived(f, 10);

      Domain x(1.0);

      Range y;
      f.evaluate(x,y);
      std::cout << "Value of f at x=(" << x <<") is " << y << std::endl;

      auto df = Dune::Functions::derivative(f);
      DRange dy;
      df->evaluate(x,dy);
      std::cout << "Value of f' at x=(" << x <<") is " << dy << std::endl;

      auto ddf = Dune::Functions::derivative(df);
      DDRange ddy;
      ddf->evaluate(x,ddy);
      std::cout << "Value of f'' at x=" << x <<" is " << ddy << std::endl;
    }
    return passed;
  }
};





int main ( int argc, char **argv )
try
{
  bool passed = true;

  passed = passed and FunctionFromCallableTest::check();

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
