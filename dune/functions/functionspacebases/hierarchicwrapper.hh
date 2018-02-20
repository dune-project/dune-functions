// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICWRAPPER_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICWRAPPER_HH

#include <dune/functions/functionspacebases/hierarchicaccess.hh>
#include <dune/functions/functionspacebases/hierarchicresizer.hh>

namespace Dune { namespace Functions
{
  /// \brief A wrapper providing multiindex access to matrix entries
  /**
   * Wraps a coefficient matrix of possibly hierarchic block structure
   * and provides access via multiindices in the `operator()`. This allows
   * to use the wrapped matrix in assembling methods where
   * multiindices are given by the local-indexset.
   *
   * The row and column blocking hierarchies are described by the template parameters `RowBlocking`
   * and `ColBlocking`, respectively, that should be either `Blocking::Flat`,
   * `Blocking::LeafBlocked` or `Blocking::Blocked` where the last one provides a blocking
   * structure for its childs. The blocking structure can be obtained from the GlobalBasis it is
   * based on, using \ref Blocking_t.
   *
   * \tparam Vector      Container with classical (block) access by an `operator[]`
   * \tparam RowBlocking Type from namespace `Dune::Functions::Blocking` describing the
   *                     hierarchic blocking structure of the rows.
   * \tparam ColBlocking Type from namespace `Dune::Functions::Blocking` describing the
   *                     hierarchic blocking structure of the columns.
   **/
  template <class Matrix, class RowBlocking, class ColBlocking>
  class HierarchicMatrixWrapper
  {
  public:
    /// \brief Constructor, stores a reference to the `matrix`.
    explicit HierarchicMatrixWrapper(Matrix& matrix)
      : matrix_(matrix)
      , resizer_(matrix)
    {}

    /// \brief Mutable access to matrix entries using multiindex access
    template <class RowIndex, class ColIndex>
    decltype(auto) operator()(RowIndex const& row, ColIndex const& col)
    {
      return matrixAccess(matrix_, RowBlocking{}, ColBlocking{}, row, col);
    }

    /// \brief Const access to matrix entries using multiindex access
    template <class RowIndex, class ColIndex>
    decltype(auto) operator()(RowIndex const& row, ColIndex const& col) const
    {
      return matrixAccess(matrix_, RowBlocking{}, ColBlocking{}, row, col);
    }

    /// \brief Resize all blocks in the hierarchy using the \ref SizeInfo size containers.
    template <class RowSize, class ColSize>
    void resize(RowSize const& rowSize, ColSize const& colSize)
    {
      resizer_(rowSize, colSize);
    }

  private:
    Matrix& matrix_;
    HierarchicMatrixResizer<Matrix, RowBlocking, ColBlocking> resizer_;
  };

  /// \brief Generator for \ref HierarchicMatrixWrapper. \relates HierarchicMatrixWrapper
  template <class Matrix, class RowBlocking, class ColBlocking>
  auto hierarchicMatrixWrapper(Matrix& matrix, RowBlocking const&, ColBlocking const&)
  {
    return HierarchicMatrixWrapper<Matrix, RowBlocking, ColBlocking>{matrix};
  }


  /// \brief A wrapper providing multiindex access to vector entries
  /**
   * Wraps a coefficient vector of possibly hierarchic block structure
   * and provides access via multiindices in the `operator[]`. This allows
   * to use the wrapped vector in interpolation and assembling methods where
   * multiindices are given by the local-indexset.
   *
   * The blocking hierarchy is described by the template parameter `BlockingTag`
   * that should be either `Blocking::Flat`, `Blocking::LeafBlocked` or `Blocking::Blocked`
   * where the last one provides a blocking structure for its childs. The blocking
   * structure can be obtained from the GlobalBasis it is based on, using \ref Blocking_t.
   *
   * \tparam Vector      Container with classical (block) access by an `operator[]`
   * \tparam BlockingTag Type from namespace `Dune::Functions::Blocking` describing the
   *                     hierarchic blocking structure.
   **/
  template <class Vector, class BlockingTag>
  class HierarchicVectorWrapper
  {
  public:
    /// \brief Constructor, stores a reference to the `vector`.
    explicit HierarchicVectorWrapper(Vector& vector)
      : vector_(vector)
      , resizer_(vector)
    {}

    /// \brief Mutable access to vector entries using multiindex access
    template <class MultiIndex>
    decltype(auto) operator[](MultiIndex const& idx)
    {
      return vectorAccess(vector_, BlockingTag{}, idx);
    }

    /// \brief Const access to vector entries using multiindex access
    template <class MultiIndex>
    decltype(auto) operator[](MultiIndex const& idx) const
    {
      return vectorAccess(vector_, BlockingTag{}, idx);
    }

    /// \brief Resize all blocks in the hierarchy using the \ref SizeInfo size container.
    template <class SizeInfo>
    void resize(SizeInfo const& sizeInfo)
    {
      resizer_(sizeInfo);
    }

  private:
    Vector& vector_;
    HierarchicVectorResizer<Vector, BlockingTag> resizer_;
  };

  /// \brief Generator for \ref HierarchicVectorWrapper. \relates HierarchicVectorWrapper
  template <class Vector, class BlockingTag>
  auto hierarchicVectorWrapper(Vector& vector, BlockingTag const&)
  {
    return HierarchicVectorWrapper<Vector, BlockingTag>{vector};
  }

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICWRAPPER_HH
