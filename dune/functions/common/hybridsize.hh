// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_HYBRIDSIZE_HH
#define DUNE_FUNCTIONS_COMMON_HYBRIDSIZE_HH

#include <utility>

#include <dune/common/indices.hh>
#include <dune/functions/common/access.hh>


namespace Dune
{
  // forward declarations
  template <class... T>
  class TupleVector;

  template <class FirstRow, class... Args>
  class MultiTypeBlockMatrix;

  template <class... Args>
  class MultiTypeBlockVector;
}

namespace Dune { namespace Functions
{
  namespace Concept
  {
    template <class V>
    static constexpr bool DynamicVectorAccessible()
    {
      return VectorAccessible<V, std::size_t>();
    }

    template <class M>
    static constexpr bool DynamicMatrixAccessible()
    {
      return MatrixAccessible<M, std::size_t, std::size_t>() ||
             Concept::isCallable<M, std::size_t, std::size_t>();
    }

  } // end namespace Concept


  namespace Impl
  {
    template <class Matrix>
    struct SizeImpl;

    template <class... Ts>
    struct SizeImpl<TupleVector<Ts...>>
        : public index_constant<sizeof...(Ts)> {};

    template <class... Ts>
    struct SizeImpl<MultiTypeBlockVector<Ts...>>
        : public index_constant<sizeof...(Ts)> {};

  } // end namespace Impl


#ifdef DOXYGEN
  /// \brief Return the size of the vector `vec`.
  /**
   * If the vector `vec` can be accessed using (dynamic) indices, the function returns the number
   * of entries as integer value. Otherwise a `std::integral_constant` is returned.
   **/
  template <class Vector>
  implementation-defined hybridSize(Vector const& vec);
#else

  template <class Vector,
    std::enable_if_t<Concept::DynamicVectorAccessible<Vector>(), int> = 0>
  auto hybridSize(Vector const& vec) -> decltype(vec.size()) { return vec.size(); }

  template <class Vector,
    std::enable_if_t<not Concept::DynamicVectorAccessible<Vector>(), int> = 0>
  constexpr typename Impl::SizeImpl<Vector>::type hybridSize(Vector const& vec) { return {}; }

#endif

  namespace Impl
  {
    template <class Matrix>
    struct NumRowsImpl;

    template <class... Rows>
    struct NumRowsImpl<MultiTypeBlockMatrix<Rows...>>
        : public index_constant<sizeof...(Rows)> {};

  } // end namespace Impl


#ifdef DOXYGEN
  /// \brief Return the number of rows of the matrix `mat`.
  /**
   * If the matrix `mat` can be accessed using (dynamic) indices, the function returns the number
   * of rows as integer value. Otherwise a `std::integral_constant` is returned.
   **/
  template <class Matrix>
  implementation-defined hybridNumRows(Matrix const& mat);
#else

  template <class Matrix,
    std::enable_if_t<Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  auto hybridNumRows(Matrix const& mat) -> decltype(mat.num_rows()) { return mat.num_rows(); }

  template <class Matrix,
    std::enable_if_t<Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  auto hybridNumRows(Matrix const& mat) -> decltype(mat.N()) { return mat.N(); }

  template <class Matrix,
    std::enable_if_t<not Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  constexpr typename Impl::NumRowsImpl<Matrix>::type hybridNumRows(Matrix const& mat) { return {}; }

#endif


  namespace Impl
  {
    template <class Matrix>
    struct NumColsImpl;

    template <class FirstRow, class... Rows>
    struct NumColsImpl<MultiTypeBlockMatrix<FirstRow,Rows...>>
        : public SizeImpl<FirstRow> {};

  } // end namespace Impl


#ifdef DOXYGEN
  /// \brief Return the number of columns of the matrix `mat`.
  /**
   * If the matrix `mat` can be accessed using (dynamic) indices, the function returns the number
   * of columns as integer value. Otherwise a `std::integral_constant` is returned.
   **/
  template <class Matrix>
  implementation-defined hybridNumCols(Matrix const& mat);
#else

  template <class Matrix,
    std::enable_if_t<Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  auto hybridNumCols(Matrix const& mat) -> decltype(mat.num_rows()) { return mat.num_cols(); }

  /// \brief Return the number of columns of a matrix as integer value
  template <class Matrix,
    std::enable_if_t<Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  auto hybridNumCols(Matrix const& mat) -> decltype(mat.M()) { return mat.M(); }

  /// \brief Return integral_constant representing the number of columns of a matrix
  template <class Matrix,
    std::enable_if_t<not Concept::DynamicMatrixAccessible<Matrix>(), int> = 0>
  constexpr typename Impl::NumColsImpl<Matrix>::type hybridNumCols(Matrix const& mat) { return {}; }

#endif

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_HYBRIDSIZE_HH
