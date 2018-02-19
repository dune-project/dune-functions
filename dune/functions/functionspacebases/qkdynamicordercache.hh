// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_QK_DYNAMIC_ORDER_CACHE_HH
#define DUNE_FUNCTIONS_QK_DYNAMIC_ORDER_CACHE_HH

#include <map>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>
#include <dune/localfunctions/lagrange/qk.hh>

#include <dune/localfunctions/lagrange/p0.hh>


/** \brief A cache that holds local finite elements of different orders (created by
 * a Factory class), dynamically
 */
namespace Dune {
namespace Functions {
  template<class D, class R, int dim, class Factory>
  class DynamicOrderLocalFiniteElementCache
  {
  protected:
    using FEFactory = Factory;
    using T = typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits;
    using FE = LocalFiniteElementVirtualInterface<T>;
    using FEMap = typename std::map<int, FE*>;

  public:
    /** \brief Type of the finite elements stored in this cache */
    using FiniteElementType = FE;

    /** \brief Default constructor */
    DynamicOrderLocalFiniteElementCache() {}

    /** \brief Copy constructor */
    DynamicOrderLocalFiniteElementCache(const DynamicOrderLocalFiniteElementCache& other)
    {
      auto it = other.cache_.begin();
      auto end = other.cache_.end();
      for(; it!=end; ++it)
        cache_[it->first] = (it->second)->clone();
    }

    ~DynamicOrderLocalFiniteElementCache()
    {
      auto it = cache_.begin();
      auto end = cache_.end();
      for(; it!=end; ++it)
        delete it->second;
    }

    //! Get local finite element for given order
    const FiniteElementType& get(const int& degree) const
    {
      auto it = cache_.find(degree);
      if (it==cache_.end())
      {
        cache_[degree] = FEFactory::create(degree);

        if (cache_[degree]==nullptr)
          DUNE_THROW(Dune::NotImplemented,"No local finite element available for order " << degree);

        return *(cache_[degree]);
      }
      return *(it->second);
    }

  protected:
    mutable FEMap cache_;

  };

  /**
   * \brief Factory for Qk finite elements (using equidistant Lagrange nodes)
   */
  template<class D, class R, int dim>
  struct DynamicOrderQkLocalFiniteElementFactory
  {
    typedef typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;

    template<int k>
    using Qk = Dune::QkLocalFiniteElement<D, R, dim, k>;


    //! create finite element for given GeometryType
    static FiniteElementType* create(const size_t& gt)
    {
      switch(gt) {
        case(0):
          return new  LocalFiniteElementVirtualImp<Qk<0>>(Qk<0>());
        case(1):
          return new  LocalFiniteElementVirtualImp<Qk<1>>(Qk<1>());
        case(2):
          return new  LocalFiniteElementVirtualImp<Qk<2>>(Qk<2>());
        case(3):
          return new  LocalFiniteElementVirtualImp<Qk<3>>(Qk<3>());
        case(4):
          return new  LocalFiniteElementVirtualImp<Qk<4>>(Qk<4>());
        case(5):
          return new  LocalFiniteElementVirtualImp<Qk<5>>(Qk<5>());
        case(6):
          return new  LocalFiniteElementVirtualImp<Qk<6>>(Qk<6>());
        case(7):
          return new  LocalFiniteElementVirtualImp<Qk<7>>(Qk<7>());
        case(8):
          return new  LocalFiniteElementVirtualImp<Qk<8>>(Qk<8>());
        case(9):
          return new  LocalFiniteElementVirtualImp<Qk<9>>(Qk<9>());
        case(10):
          return new  LocalFiniteElementVirtualImp<Qk<10>>(Qk<10>());
        case(11):
          return new  LocalFiniteElementVirtualImp<Qk<11>>(Qk<11>());
        case(12):
          return new  LocalFiniteElementVirtualImp<Qk<12>>(Qk<12>());
        case(13):
          return new  LocalFiniteElementVirtualImp<Qk<13>>(Qk<13>());
        case(14):
          return new  LocalFiniteElementVirtualImp<Qk<14>>(Qk<14>());
        default:
          DUNE_THROW(Dune::NotImplemented, "Dynamic Qk only up to degree 14");
      }
    }
  };

  /**
   * \brief Alias for a dynamic order Qk FE cache.
   */
  template<class D, class R, int dim>
  using DynamicQkFiniteElementCache = DynamicOrderLocalFiniteElementCache<D, R, dim, DynamicOrderQkLocalFiniteElementFactory<D, R, dim>>;
}
}
#endif
