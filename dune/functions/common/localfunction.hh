// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_LOCALFUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_LOCALFUNCTION_HH

#include <memory>
#include <dune/functions/common/differentiablefunction.hh>

namespace Dune {

namespace Functions {

  /**
   * \tparam DF Type derived from DifferentiableFunction
   * \tparam LC Type of local context, e.g. grid cell
   * \tparam R  Range type
   */
template<typename DF, typename LC>
class LocalFunction
  : public DifferentiableFunction<typename LC::Geometry::LocalCoordinate, typename DF::Range>
{

  typedef DifferentiableFunction<typename LC::Geometry::LocalCoordinate, typename DF::Range> Base;

public:

  typedef DF GlobalFunction;
  typedef LC LocalContext;

  typedef typename LC::Geometry::LocalCoordinate Domain;

  typedef LocalFunction<typename DF::Derivative,LC> Derivative;
  typedef std::shared_ptr<typename Base::Derivative> DerivativeBasePointer;

  virtual void bind(const LocalContext&) = 0;

  virtual void unbind() = 0;

  virtual DerivativeBasePointer derivative() const = 0;

  virtual const LocalContext& localContext() const = 0;

};


template<typename Function>
shared_ptr<typename Function::LocalFunction> localFunction(const Function& f)
{
  return std::static_pointer_cast<typename Function::LocalFunction>(f.localFunction());
}

template<typename Function>
shared_ptr<typename Function::LocalFunction> localFunction(const Function& f, const typename Function::Element& e)
{
  auto p = std::static_pointer_cast<typename Function::LocalFunction>(f.localFunction());
  p->bind(e);
  return std::move(p);
}


} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_LOCALFUNCTION_HH
