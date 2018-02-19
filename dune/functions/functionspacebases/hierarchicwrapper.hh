#pragma once

#include <dune/functions/functionspacebases/hierarchicaccess.hh>
#include <dune/functions/functionspacebases/hierarchicresizer.hh>

namespace Dune { namespace Functions
{
  template <class Matrix, class RowBlocking, class ColBlocking>
  class HierarchicMatrixWrapper
  {
  public:
    explicit HierarchicMatrixWrapper(Matrix& matrix)
      : matrix_(matrix)
      , resizer_(matrix)
    {}

    template <class RowIndex, class ColIndex>
    decltype(auto) operator()(RowIndex const& row, ColIndex const& col)
    {
      return matrixAccess(matrix_, RowBlocking{}, ColBlocking{}, row, col);
    }

    template <class RowSize, class ColSize>
    void resize(RowSize const& rowSize, ColSize const& colSize)
    {
      resizer_(rowSize, colSize);
    }

  private:
    Matrix& matrix_;
    HierarchicMatrixResizer<Matrix, RowBlocking, ColBlocking> resizer_;
  };

  template <class Matrix, class RowBlocking, class ColBlocking>
  auto hierarchicMatrixWrapper(Matrix& matrix, RowBlocking const&, ColBlocking const&)
  {
    return HierarchicMatrixWrapper<Matrix, RowBlocking, ColBlocking>{matrix};
  }


  template <class Vector, class BlockingTag>
  class HierarchicVectorWrapper
  {
  public:
    explicit HierarchicVectorWrapper(Vector& vector)
      : vector_(vector)
      , resizer_(vector)
    {}

    template <class Index>
    decltype(auto) operator[](Index const& idx)
    {
      return vectorAccess(vector_, BlockingTag{}, idx);
    }

    template <class Index>
    decltype(auto) operator[](Index const& idx) const
    {
      return vectorAccess(vector_, BlockingTag{}, idx);
    }

    template <class SizeInfo>
    void resize(SizeInfo const& sizeInfo)
    {
      resizer_(sizeInfo);
    }

  private:
    Vector& vector_;
    HierarchicVectorResizer<Vector, BlockingTag> resizer_;
  };

  template <class Vector, class BlockingTag>
  auto hierarchicVectorWrapper(Vector& vector, BlockingTag const&)
  {
    return HierarchicVectorWrapper<Vector, BlockingTag>{vector};
  }

}} // end namespace Dune::Functions
