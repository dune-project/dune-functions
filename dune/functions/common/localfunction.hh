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
   * \tparam LC  Type of local context, e.g. grid cell
   */
template<typename DF, typename LC>
class LocalFunction
  : public DifferentiableFunction<typename LC::Geometry::LocalCoordinate, typename DF::Range>
{

  typedef DF GlobalFunction;
  typedef LC LocalContext;

  typedef typename LC::Geometry::LocalCoordinate Domain;

  typedef LocalFunction<Domain,typename DF::DerivativeRange> Derivative;

  virtual void bind(const LocalContext&) = 0;

  virtual void unbind() = 0;

  virtual Derivative* derivative() = 0;

  virtual const LocalContext& localContext() = 0;

};


template<typename Function>
shared_ptr<typename Function::ElementFunction> elementFunction(const Function& f)
{
  return static_pointer_cast<typename Function::ElementFunction>(f.elementFunction());
}

template<typename Function>
shared_ptr<typename Function::ElementFunction> elementFunction(const Function& f, const typename Function::Element& e)
{
  auto p = static_pointer_cast<typename Function::ElementFunction>(f.elementFunction());
  p.bind(e);
  return std::move(p);
}


} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_LOCALFUNCTION_HH
