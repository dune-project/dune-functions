// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TEST_DERIVATIVE_CHECK_HH
#define DUNE_FUNCTIONS_COMMON_TEST_DERIVATIVE_CHECK_HH

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/functions/common/differentiablefunction.hh>


template<class F>
struct DerivativeCheck
{

  static bool checkType(const F& f)
  {
    bool passed = true;

    typedef typename Dune::Functions::DerivativeTraits<typename F::Domain, typename F::Range>::DerivativeRange DR;
    typedef typename Dune::Functions::DifferentiableFunction<typename F::Domain, DR> DerivativeInterface;

    // DerivativeRange should match DerivativeTraits::DerivativeRange
    passed = passed and Dune::is_same<DR, typename F::DerivativeRange>::value;

    // Derivative should be derived from interface
    passed = passed and Dune::IsBaseOf<DerivativeInterface, typename F::Derivative>::value;

    return passed;
  }

  static bool isTrulyDerived(const F& f)
  {
    typedef typename Dune::Functions::DerivativeTraits<typename F::Domain, typename F::Range>::DerivativeRange DR;
    typedef typename Dune::Functions::DifferentiableFunction<typename F::Domain, DR> DerivativeInterface;

    // To benefit from 'final' all implementation in dune-functions should export the derived type.
    // For usercode it's ok to use the interface if you willing to pay the price
    return not(Dune::is_same<DerivativeInterface, typename F::Derivative>::value);
  }

  static bool isImplemented(const F& f)
  {
    try
    {
      auto df = Dune::Functions::derivative(f);
    }
    catch (Dune::NotImplemented)
    {
      return false;
    }
    return true;
  }

  static bool checkAllImplementatedTrulyDerived(const F& f, int maxOrder, int order=1)
  {
    if (order>maxOrder)
    {
      std::cout << "Stopping derivative check after maxOrder=" << maxOrder << std::endl;
      return true;
    }
    std::cout << "Checking derivative of order=" << order << std::endl;
    bool implemented = isImplemented(f);
    bool derived = isTrulyDerived(f);

    if (implemented and derived)
    {
      auto df = Dune::Functions::derivative(f).shared_ptr();
      return checkType(f) and DerivativeCheck<typename F::Derivative>::checkAllImplementatedTrulyDerived(*df, maxOrder, order+1);
    }
    if (implemented and not(derived))
    {
      std::cout << "Derivative of order=" << order << " is implemented but typedef ::Derivative is the interface class." << std::endl;
      return false;
    }
    if (not(implemented) and derived)
    {
      std::cout << "Warning: Typedef ::Derivative of order=" << order << " is a derived type but derivative is not implemented." << std::endl;
    }
    return true;
  }

};

#endif // DUNE_FUNCTIONS_COMMON_TEST_DERIVATIVE_CHECK_HH
