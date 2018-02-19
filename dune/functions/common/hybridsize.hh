#pragma once

#include <utility>

#include <dune/common/indices.hh>
#include <dune/common/tuplevector.hh>
#include <dune/functions/common/access.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#endif

namespace Dune { namespace Functions
{
  namespace Concept
  {
    template <class V>
    static constexpr bool DynamicVectorAccessible() { return VectorAccessible<V, std::size_t>(); }

    template <class M>
    static constexpr bool DynamicMatrixAccessible() { return MatrixAccessible<M, std::size_t, std::size_t>() || Concept::isCallable<M, std::size_t, std::size_t>(); }

  } // end namespace Concept


  namespace Impl
  {
    template <class Matrix>
    struct SizeImpl;

    template <class... Ts>
    struct SizeImpl<TupleVector<Ts...>>
        : public index_constant<sizeof...(Ts)> {};

#if HAVE_DUNE_ISTL
    template <class... Ts>
    struct SizeImpl<MultiTypeBlockVector<Ts...>>
        : public index_constant<sizeof...(Ts)> {};
#endif
  }

  // access i'th component of a vector
  template <class Vector,
    std::enable_if_t<Concept::DynamicVectorAccessible<Vector>(), int> = 0>
  auto hybridSize(Vector const& vec) { return vec.size(); }

  // return integral_constant representing the size
  template <class Vector,
    std::enable_if_t<not Concept::DynamicVectorAccessible<Vector>(), int> = 0>
  constexpr typename Impl::SizeImpl<Vector>::type hybridSize(Vector const& vec) { return {}; }

//   template <class Vector,
//     std::enable_if_t<not Concept::DynamicVectorAccessible<Vector>(), int> = 0>
//   constexpr auto hybridSize(Vector const& vec) { return Hybrid::size(vec); }


  namespace Impl
  {
    template <class Matrix>
    struct NumRowsImpl;

#if HAVE_DUNE_ISTL
    template <class... Rows>
    struct NumRowsImpl<MultiTypeBlockMatrix<Rows...>>
        : public index_constant<sizeof...(Rows)> {};
#endif
  }

  // return the number of rows of a matrix
  template <class Matrix,
    std::enable_if_t<Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  auto hybridNumRows(Matrix const& mat) -> decltype(mat.num_rows()) { return mat.num_rows(); }

  template <class Matrix,
    std::enable_if_t<Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  auto hybridNumRows(Matrix const& mat) -> decltype(mat.N()) { return mat.N(); }

  // return integral_constant representing the number of rows
  template <class Matrix,
    std::enable_if_t<not Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  constexpr typename Impl::NumRowsImpl<Matrix>::type hybridNumRows(Matrix const& mat) { return {}; }


  namespace Impl
  {
    template <class Matrix>
    struct NumColsImpl;

#if HAVE_DUNE_ISTL
    template <class FirstRow, class... Rows>
    struct NumColsImpl<MultiTypeBlockMatrix<FirstRow,Rows...>>
        : public SizeImpl<FirstRow> {};
#endif
  }

  // return the number of columns of a matrix
  template <class Matrix,
    std::enable_if_t<Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  auto hybridNumCols(Matrix const& mat) -> decltype(mat.num_rows()) { return mat.num_cols(); }

  template <class Matrix,
    std::enable_if_t<Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  auto hybridNumCols(Matrix const& mat) -> decltype(mat.M()) { return mat.M(); }

  // return integral_constant representing the number of columns
  template <class Matrix,
    std::enable_if_t<not Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  constexpr typename Impl::NumColsImpl<Matrix>::type hybridNumCols(Matrix const& mat) { return {}; }

}} // end namespace Dune::Functions
