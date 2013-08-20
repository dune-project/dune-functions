// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DISCRETEGRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_DISCRETEGRIDVIEWFUNCTION_HH

#include <memory>
#include <dune/functions/common/gridviewfunction.hh>

namespace Dune {

namespace Functions {

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

#endif // DUNE_FUNCTIONS_COMMON_DISCRETEGRIDVIEWFUNCTION_HH
