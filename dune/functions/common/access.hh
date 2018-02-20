// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_ACCESS_HH
#define DUNE_FUNCTIONS_COMMON_ACCESS_HH

#include <utility>
#include <dune/functions/common/functionconcepts.hh>

namespace Dune { namespace Functions
{
  namespace Concept
  {
    struct HasVectorAccess
    {
      template <class V, class I>
      auto require(V&& v, I&& i) -> decltype( v[i] );
    };

    struct HasMatrixAccess
    {
      template <class M, class I, class J>
      auto require(M&& m, I&& i, J&& j) -> decltype( m[i][j] );
    };

    template <class V, class I>
    static constexpr bool VectorAccessible() { return models<HasVectorAccess, V, I>(); }

    template <class M, class I, class J>
    static constexpr bool MatrixAccessible() { return models<HasMatrixAccess, M, I, J>(); }

  } // end namespace Concept


#ifdef DOXYGEN
  /// \brief Uniform vector access using [.]
  template <class Vector, class I>
  decltype(auto) access(Vector&& vec, I const& i);

  /// \brief Uniform matrix access using either [.][.] or (.,.)
  template <class Matrix, class I, class J>
  decltype(auto) access(Matrix&& mat, I const& i, J const& j);
#else

  // access i'th component of a vector
  template <class Vector, class I,
    std::enable_if_t<Concept::VectorAccessible<Vector,I>(), int> = 0>
  decltype(auto) access(Vector&& vec, I const& i) { return vec[i]; }

  // fall-back implementation for scalars
  template <class Vector, class I,
    std::enable_if_t<not Concept::VectorAccessible<Vector,I>(), int> = 0>
  decltype(auto) access(Vector&& vec, I const& /*i*/) { return std::forward<Vector>(vec); }


  // access (i,j)'th component of a matrix using [.][.]
  template <class Matrix, class I, class J,
    std::enable_if_t<Concept::MatrixAccessible<Matrix,I,J>(), int> = 0>
  decltype(auto) access(Matrix&& mat, I const& i, J const& j) { return mat[i][j]; }

  // access (i,j)'th component of a matrix using (.,.)
  template <class Matrix, class I, class J,
    std::enable_if_t<not Concept::MatrixAccessible<Matrix,I,J>() && Concept::isCallable<Matrix,I,J>(), int> = 0>
  decltype(auto) access(Matrix&& mat, I const& i, J const& j) { return mat(i,j); }

  // fall-back implementation for scalars
  template <class Matrix, class I, class J,
    std::enable_if_t<not Concept::MatrixAccessible<Matrix,I,J>() && not Concept::isCallable<Matrix,I,J>(), int> = 0>
  decltype(auto) access(Matrix&& mat, I const& /*i*/, J const& /*j*/) { return std::forward<Matrix>(mat); }

#endif

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_ACCESS_HH
