// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MATRIXVECTORGENERATOR_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MATRIXVECTORGENERATOR_HH

#include <memory>
#include <utility>
#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/functions/functionspacebases/blocking.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#else
#error "Need dune-istl for this implementation"
#endif

namespace Dune
{
  // forward declarations
  template <class FirstRow, class... Args>
  class MultiTypeBlockMatrix;

  template <class... Args>
  class MultiTypeBlockVector;
}


namespace Dune { namespace Functions
{
  namespace Impl
  {
    template <class Value, template <class...> class Container, class BlockingTag, class = void>
    struct VectorGenerator;

    template <class Value, template <class...> class Container, class RowBlocking, class ColBlocking, class = void>
    struct MatrixGenerator;

  } // end namespace Impl

  /// \brief Generate a vector container compatible with the blocking structure
  template <class Value, template <class...> class Container, class BlockingTag>
  using VectorGenerator_t = typename Impl::VectorGenerator<Value, Container, BlockingTag>::type;

  /// \brief Generate a matrix container compatible with the blocking structure
  template <class Value, template <class...> class Container, class RowBlocking, class ColBlocking>
  using MatrixGenerator_t = typename Impl::MatrixGenerator<Value, Container, RowBlocking, ColBlocking>::type;


  namespace Impl
  {
    template <class... Ts>
    struct Types {};


    template <class... Ts> struct IsSame;

    template <class T0, class... Ts>
    struct IsSame<T0, Ts...>
        : public std::is_same<Types<T0,Ts...>, Types<Ts...,T0>> {};

    template <>
    struct IsSame<> { enum { value = true }; };


    // implementation details of the vector container
    template <class Value, template <class...> class Container, class... T>
    struct VectorGenerator<Value, Container, Blocking::Blocked<T...>, std::enable_if_t<not IsSame<T...>::value> >
    {
      using type = MultiTypeBlockVector<VectorGenerator_t<Value, Container, T>...>;
    };

    template <class Value, template <class...> class Container, class T0, class... T>
    struct VectorGenerator<Value, Container, Blocking::Blocked<T0,T...>, std::enable_if_t<IsSame<T0,T...>::value> >
    {
      using type = BlockVector<VectorGenerator_t<Value, Container, T0>>;
    };

    template <class Value, template <class...> class Container>
    struct VectorGenerator<Value, Container, Blocking::Flat>
    {
      using type = Container<Value>;
    };

    template <class Value, template <class...> class Container, std::size_t N>
    struct VectorGenerator<Value, Container, Blocking::LeafBlocked<N>>
    {
      using type = Container<FieldVector<Value,int(N)>>;
    };


    // implementation details for the matrix-container
    template <class Value, template <class...> class Container, class... S, class... T>
    struct MatrixGenerator<Value, Container, Blocking::Blocked<S...>, Blocking::Blocked<T...>,
      std::enable_if_t<not(IsSame<S...>::value && IsSame<T...>::value)> >
    {
      template <class Si>
      using Row = MultiTypeBlockVector<MatrixGenerator_t<Value, Container, Si, T>...>;

      using type = MultiTypeBlockMatrix<Row<S>...>;
    };

    template <class Value, template <class...> class Container, class S0, class... S, class T0, class... T>
    struct MatrixGenerator<Value, Container, Blocking::Blocked<S0,S...>, Blocking::Blocked<T0,T...>,
      std::enable_if_t<IsSame<S0,S...>::value && IsSame<T0,T...>::value> >
    {
      using type = Matrix<MatrixGenerator_t<Value, Container, S0, T0>>;
    };

    template <class Value, template <class...> class Container>
    struct MatrixGenerator<Value, Container, Blocking::Flat, Blocking::Flat>
    {
      using type = Container<Value>;
    };

    template <class Value, template <class...> class Container, std::size_t N, std::size_t M>
    struct MatrixGenerator<Value, Container, Blocking::LeafBlocked<N>, Blocking::LeafBlocked<M>>
    {
      using type = Container<FieldMatrix<Value,int(N),int(M)>>;
    };

    template <class Value, template <class...> class Container, std::size_t N>
    struct MatrixGenerator<Value, Container, Blocking::LeafBlocked<N>, Blocking::Flat>
    {
      using type = Container<FieldMatrix<Value,int(N),1>>;
    };

    template <class Value, template <class...> class Container, std::size_t M>
    struct MatrixGenerator<Value, Container, Blocking::Flat, Blocking::LeafBlocked<M>>
    {
      using type = Container<FieldMatrix<Value,1,int(M)>>;
    };

    template <class Value, template <class...> class Container, class... T>
    struct MatrixGenerator<Value, Container, Blocking::Blocked<T...>, Blocking::Flat, std::enable_if_t<not IsSame<T...>::value> >
    {
      using type = MultiTypeBlockMatrix<MultiTypeBlockVector<MatrixGenerator_t<Value, Container, T, Blocking::Flat>>...>;
    };

    template <class Value, template <class...> class Container, class T0, class... T>
    struct MatrixGenerator<Value, Container, Blocking::Blocked<T0,T...>, Blocking::Flat, std::enable_if_t<IsSame<T0,T...>::value> >
    {
      using type = Matrix<MatrixGenerator_t<Value, Container, T0, Blocking::Flat>>;
    };

    template <class Value, template <class...> class Container, class... T>
    struct MatrixGenerator<Value, Container, Blocking::Flat, Blocking::Blocked<T...>, std::enable_if_t<not IsSame<T...>::value> >
    {
      using type = MultiTypeBlockMatrix<MultiTypeBlockVector<MatrixGenerator_t<Value, Container, Blocking::Flat, T>...>>;
    };

    template <class Value, template <class...> class Container, class T0, class... T>
    struct MatrixGenerator<Value, Container, Blocking::Flat, Blocking::Blocked<T0,T...>, std::enable_if_t<IsSame<T0,T...>::value> >
    {
      using type = Matrix<MatrixGenerator_t<Value, Container, Blocking::Flat, T0>>;
    };

  } // end namespace Impl

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MATRIXVECTORGENERATOR_HH
