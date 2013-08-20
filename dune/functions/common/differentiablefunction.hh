// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
#include <memory>
#include <dune/common/function.hh>

namespace Dune {

namespace Functions {



template<class DT, class RT> class DifferentiableFunction;

template<class FImp> class FunctionHandle;



class InvalidRange
{};


template<class DT, class RT>
struct DerivativeTraits
{
    typedef InvalidRange DerivativeRange;
};

template<>
struct DerivativeTraits<double, double>
{
    typedef double DerivativeRange;
};



template<class DT, class RT>
class DifferentiableFunction :
    public Dune::VirtualFunction<DT, RT>,
    public std::enable_shared_from_this<DifferentiableFunction<DT, RT> >
{
    public:
        typedef DT Domain;
        typedef RT Range;
        typedef typename DerivativeTraits<DT, RT>::DerivativeRange DerivativeRange;

        typedef DifferentiableFunction<Domain, DerivativeRange> Derivative;

    private:
  //template<typename T>
  //    friend FunctionHandle<typename T::Derivative> ::Dune::Functions::derivative(const T&);

    protected:
    public:

        virtual Derivative* derivative() const = 0;

};


template<class FImp>
class FunctionHandle :
    public DifferentiableFunction<typename FImp::Domain, typename FImp::Range>
{
        typedef DifferentiableFunction<typename FImp::Domain, typename FImp::Range> Base;
    public:
        typedef typename Base::Range Range;
        typedef typename Base::Domain Domain;
        typedef typename Base::DerivativeRange DerivativeRange;

        typedef typename FImp::Derivative Derivative;

        // to be dicussed
        typedef FImp HandledFunction;
        typedef typename std::shared_ptr<const FImp> SharedPtr;

        explicit FunctionHandle(const FImp* f) :
            f_(f)
        {}


        virtual void evaluate(const Domain& x, Range& y) const DUNE_FINAL
        {
            f_->evaluate(x, y);
        }

        std::shared_ptr<const FImp> shared_ptr() const
        {
            return f_->shared_from_this();
        }

    private:
  // template<typename T>
  //    friend FunctionHandle<typename T::Derivative> ::Dune::Functions::derivative(const T&);

    protected:
    public:

        virtual Derivative* derivative() const DUNE_FINAL
        {
            return f_->derivative();
        }

        const FImp* f_;
};


template<class FImp>
FunctionHandle<typename FImp::Derivative> derivative(const FImp& f)
{
    return FunctionHandle<typename FImp::Derivative>(f.derivative());
}






template<class DT, class RT>
class InvalidFunction :
    Dune::VirtualFunction<DT, RT>
{
        typedef typename Dune::VirtualFunction<DT, RT> Base;
    public:
        typedef typename Base::Range Range;
        typedef typename Base::Domain Domain;
        typedef typename Base::DerivativeRange DerivativeRange;
        typedef InvalidFunction<DT, DerivativeRange> Derivative;

        virtual void evaluate(const Domain& x, Range& y) const
        {
            throw 1;
        }

    protected:
        virtual Derivative* derivative() const
        {
            throw 1;
        }

};



/**
 * derivative(f).evaluate(x,y)
 * auto df = derivative(f);
 * auto dfp = derivative(f).shared_ptr();
 */

} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
