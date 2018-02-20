// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_RESIZE_HH
#define DUNE_FUNCTIONS_COMMON_RESIZE_HH

#include <utility>
#include <type_traits>

#include <dune/common/typetraits.hh>

namespace Dune { namespace Functions
{
  namespace Impl
  {
    template <class Vector, class = void>
    struct VectorResizeImpl {};

    template <class Vector>
    struct VectorResizeImpl<Vector, void_t<decltype(std::declval<Vector>().resize(std::size_t(0)))>>
    {
      using type = void;
      static void apply(Vector& vec, std::size_t s) { vec.resize(s); }
    };

    template <class Vector>
    struct VectorResizeImpl<Vector, void_t<decltype(std::declval<Vector>().change_dim(std::size_t(0)))>>
    {
      using type = void;
      static void apply(Vector& vec, std::size_t s) { vec.change_dim(s); }
    };

  } // end namespace Impl

  namespace Concept
  {
    template <class Vector, class = void>
    struct VectorResizableImpl
        : public std::false_type {};

    template <class Vector>
    struct VectorResizableImpl<Vector, typename Impl::VectorResizeImpl<Vector>::type>
        : public std::true_type {};

    template <class Vector>
    static constexpr bool VectorResizable() { return VectorResizableImpl<Vector>::value; }

  } // end namespace Concept


#ifdef DOXYGEN
  /// \brief Uniform vector resize, using either vector.resize() or vector.change_dim()
  template <class Vector>
  void resize(Vector& vec, std::size_t size);
#else
  template <class Vector, class = typename Impl::VectorResizeImpl<Vector>::type>
  void resize(Vector& vec, std::size_t size) { Impl::VectorResizeImpl<Vector>::apply(vec,size); }
#endif


  namespace Impl
  {
    template <class Matrix, class = void>
    struct MatrixResizeImpl {};

    template <class Matrix>
    struct MatrixResizeImpl<Matrix, void_t<decltype(std::declval<Matrix>().resize(std::size_t(0),std::size_t(0)))>>
    {
      using type = void;
      static void apply(Matrix& mat, std::size_t r, std::size_t c) { mat.resize(r,c); }
    };

    template <class Matrix>
    struct MatrixResizeImpl<Matrix, void_t<decltype(std::declval<Matrix>().setSize(std::size_t(0),std::size_t(0)))>>
    {
      using type = void;
      static void apply(Matrix& mat, std::size_t r, std::size_t c) { mat.setSize(r,c); }
    };

    template <class Matrix>
    struct MatrixResizeImpl<Matrix, void_t<decltype(std::declval<Matrix>().change_dim(std::size_t(0),std::size_t(0)))>>
    {
      using type = void;
      static void apply(Matrix& mat, std::size_t r, std::size_t c) { mat.change_dim(r,c); }
    };

  } // end namespace Impl

  namespace Concept
  {
    template <class Matrix, class = void>
    struct MatrixResizableImpl
        : public std::false_type {};

    template <class Matrix>
    struct MatrixResizableImpl<Matrix, typename Impl::MatrixResizeImpl<Matrix>::type>
        : public std::true_type {};

    template <class Matrix>
    static constexpr bool MatrixResizable() { return MatrixResizableImpl<Matrix>::value; }

  } // end namespace Concept


#ifdef DOXYGEN
  /// \brief Uniform matrix resize, using either matrix.resize() or matrix.change_dim()
  template <class Matrix>
  void resize(Matrix& mat, std::size_t r, std::size_t c);
#else
  template <class Matrix, class = typename Impl::MatrixResizeImpl<Matrix>::type>
  void resize(Matrix& mat, std::size_t r, std::size_t c) { Impl::MatrixResizeImpl<Matrix>::apply(mat,r,c); }
#endif

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_RESIZE_HH
