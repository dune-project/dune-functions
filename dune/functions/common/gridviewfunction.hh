// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH

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



template<typename GV, typename RT>
class GridViewFunction
  : public DifferentiableFunction<typename GV::template Codim<0>::Geometry::GlobalCoordinate,RT>
{

  typedef DifferentiableFunction<
    typename GV::template Codim<0>::Geometry::GlobalCoordinate,
    RT
    > Base;

  typedef GV GridView;
  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;
  typedef typename Base::DerivativeRange DerivativeRange;

  typedef GridViewFunction<GV,DerivativeRange> Derivative;

  typedef typename GV::template Codim<0>::Geometry::LocalCoordinate LocalDomain;
  typedef typename GV::template Codim<0>::Entity Element;

  typedef LocalFunction<GridViewFunction,Element> ElementFunction;

protected:

  typedef shared_ptr<ElementFunction> ElementFunctionBasePointer;

public:

  virtual ElementFunctionBasePointer elementFunction() const = 0;

  // virtual void evaluate(const Element& e, const LocalDomain& coord, Range& r) const

  virtual Derivative* derivative() const = 0;

  const GridView& gridView()
  {
    return gv_;
  }

  GridViewFunction(const GridView& gv)
    : gv_(gv)
  {}

private:

  GridView gv_;

};



template<typename GFS, typename V>
class DiscreteGridViewFunction
  : public GridViewFunction<typename GFS::GlobalCoordinate,
                            typename GFS::Range
                            >
{

  class ElementFunction
    : public typename Base::ElementFunction
  {

    typedef typename Base::ElementFunction Base;

    typedef typename Base::Domain Domain;

    ElementFunction(shared_ptr<const DiscreteGridViewFunction> dgvf)
      : dgvf_(dgvf)
      , element_(nullptr)
    {}

    virtual void bind(const Element& element) DUNE_FINAL
    {
      element_ = &element;
      lfs.bind(e);
      lfs_cache.update();
      x_view.bind(lfs_cache);
      x_view.read(xl);
      x_view.unbind();
    }

    virtual void unbind() DUNE_FINAL
    {
      element_ = nullptr;
    }

    virtual void evaluate(const Domain& coord, Range& r) DUNE_FINAL
    {
      lfs.finiteElement().basis().evaluateFunction(coord,yb);
      r = 0;
      for (unsigned int i=0; i<yb.size(); i++)
        {
          r.axpy(xl[i],yb[i]);
        }
    }

    virtual Element& element() const DUNE_FINAL
    {
#ifndef NDEBUG
      if (!element_)
        DUNE_THROW(InvalidStateException,"bla");
#endif
      return *element_;
    }

  };

  virtual typename Base::ElementFunctionBasePointer elementFunction() const
  {
    return make_shared<ElementFunction>(this->make_shared_from_this());
  }

};








} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH
