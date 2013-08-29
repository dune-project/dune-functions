// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
#include <memory>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/function.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/functions/common/final.hh>

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

template<typename F, int n>
struct DerivativeTraits<FieldVector<F,n>,F>
{
  typedef FieldVector<F,n> DerivativeRange;
};

template<typename F, int n, int m>
struct DerivativeTraits<FieldVector<F,n>,FieldVector<F,m> >
{
  typedef FieldMatrix<F,n,m> DerivativeRange;
};

template<typename F, int n, int m>
struct DerivativeTraits<FieldVector<F,n>,FieldMatrix<F,m,1> >
{
  typedef FieldMatrix<F,n,m> DerivativeRange;
};



/** \brief Abstract base class for functions that allow to compute a derivative
 *
 * \tparam DT Type used for the independent variable (may be a vector type for functions of several variables)
 * \tparam RT Type used for function values
 *
 * We view differentiation as a map between functions.  In other words, the derivative of a function
 * is again a function.  The derivative function implements the same DifferentiableFunction base class,
 * only the type for function values will have changed.
 */
template<class DT, class RT>
class DifferentiableFunction :
    public Dune::VirtualFunction<DT, RT>,
    public std::enable_shared_from_this<DifferentiableFunction<DT, RT> >
{
    public:
        /** \brief Type used for the independent variable (may be a vector type for functions of several variables) */
        typedef DT Domain;

        /** \brief Type used for function values */
        typedef RT Range;

        /** \brief Type of the values of the derivative function */
        typedef typename DerivativeTraits<DT, RT>::DerivativeRange DerivativeRange;

        /** \brief Type of the derivative function */
        typedef DifferentiableFunction<Domain, DerivativeRange> Derivative;

        /** \brief Get the derivative function
         *
         * ### Life time of the return value
         *
         * This method returns a simple C pointer to the derivative object.
         * When using this pointer be aware of the following life time guarantees:
         * A priori, each function holds ownership of its derivative function object.
         * Therefore, when the function goes out of scope (or gets destroyed in some other way),
         * then its derivative (and recursively the higher-order derivatives as well) gets
         * destroyed as well.  However, DifferentiableFunction provides infrastructure for
         * being reference counted (that is why it derives from std::enable_shared_from_this).
         * In particular, on the return type of the function you can call the method
         * shared_from_this() to obtain a std::shared_ptr<Derivative>.  This std::shared_ptr then
         * shares the ownership of the derivative with the function itself.  That
         * way, the life time of a derivative can be extended beyond the life time of the
         * function it was derived from.
         *
         * If you are afraid of handling the C pointer yourself you can get the derivative object
         * in a different way.  Call the free function
         * template<class FImp> FunctionHandle<typename FImp::Derivative> derivative(const FImp& f)
         * with f being the function you want to derive.  The return value FunctionHandle
         * acts like a light-weight safe-pointer, and again allows you to obtain a std::shared_ptr
         * to the derivative.
         */
        virtual Derivative* derivative() const = 0;


        /** \brief Get a shared_ptr to this object to obtain shared ownership.
         *
         * \warning The semantics of this method are slightly different from the default version
         *          in the C++ standard. Notably, the standard version of shared_from_this() throws
         *          an exception of type std::bad_weak_ptr if no shared_ptr has been attached to the
         *          object. In contrast, this implementation will switch to non-managed mode if you
         *          call shared_from_this() without an attached shared_ptr. This change makes it
         *          possible to call shared_from_this() on stack-allocated objects.
         */
        shared_ptr<DifferentiableFunction> shared_from_this()
        {
          try
            {
              return std::enable_shared_from_this<DifferentiableFunction<DT, RT> >::shared_from_this();
            }
          catch (std::bad_weak_ptr&)
            {
              _local_ptr = stackobject_to_shared_ptr(*this);
              return _local_ptr;
            }
        }

        /** \brief Get a shared_ptr to this object to obtain shared ownership.
         *
         * \warning The semantics of this method are slightly different from the default version
         *          in the C++ standard. Notably, the standard version of shared_from_this() throws
         *          an exception of type std::bad_weak_ptr if no shared_ptr has been attached to the
         *          object. In contrast, this implementation will switch to non-managed mode if you
         *          call shared_from_this() without an attached shared_ptr. This change makes it
         *          possible to call shared_from_this() on stack-allocated objects.
         */
        shared_ptr<const DifferentiableFunction> shared_from_this() const
        {
          try
            {
              return std::enable_shared_from_this<DifferentiableFunction<DT, RT> >::shared_from_this();
            }
          catch (std::bad_weak_ptr&)
            {
              _local_ptr = stackobject_to_shared_ptr(*const_cast<DifferentiableFunction*>(this));
              return _local_ptr;
            }
        }


private:

  mutable shared_ptr<DifferentiableFunction> _local_ptr;
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

        /**
         * \brief Construct FunctionHandle from pointer to function.
         *
         * The FunctionHandle does not take ownership of this
         * pointer. So its life time has to be controlled
         * outside.
         */
        explicit FunctionHandle(const FImp* f) :
            f_(f)
        {}


        virtual void evaluate(const Domain& x, Range& y) const DUNE_FINAL
        {
            f_->evaluate(x, y);
        }

        std::shared_ptr<const FImp> shared_ptr() const
        {
            return std::static_pointer_cast<const FImp>(f_->shared_from_this());
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

/** \brief Get the derivative of a DifferentiableFunction
 *
 * \tparam FImp A DifferentiableFunction, i.e., a function that has a derivative() method
 *
 * Using this methods avoids ever having to handle the C-style pointer that is returned
 * by DifferentiableFunction::derivative.
 */
template<class FImp>
FunctionHandle<typename FImp::Derivative> derivative(const FImp& f)
{
    return FunctionHandle<typename FImp::Derivative>(f.derivative());
}


} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLEFUNCTION_HH
