#pragma once

#include <dune/common/hybridutilities.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/functions/common/hybridsize.hh>
#include <dune/functions/common/resize.hh>
#include <dune/functions/functionspacebases/blocking.hh>

namespace Dune { namespace Functions
{
  namespace Impl
  {
    template <class V, class BlockingTag, class SizeInfo, class Prefix>
    static void resizeVector(V& vec, BlockingTag b, SizeInfo const& sizeInfo, Prefix prefix)
    {
      auto s = sizeInfo.size(prefix);
      using VectorResizable = Std::bool_constant<Concept::VectorResizable<V>()>;
      using Blocked = Std::bool_constant<Concept::Blocked<BlockingTag>()>;

      Hybrid::ifElse(VectorResizable{},
        [&](auto id) { resize(id(vec), s); },
        [&](auto id) { assert(hybridSize(vec) == s && "Wrong (static) dimension."); });

      Hybrid::ifElse(Blocked{},
        [&](auto id) {
          prefix.push_back(0);
          Hybrid::forEach(range(hybridSize(vec)), [&](auto const _i) {
            prefix.back() = _i;
            resizeVector(id(vec)[_i], b[_i], sizeInfo, prefix);
          });
        });
    }

  } // end namespace Impl

  template <class Vector, class BlockingTag>
  class HierarchicVectorResizer
  {
  public:
    HierarchicVectorResizer(Vector& vector)
      : vector_(vector)
    {}

    template <class SizeInfo, class Prefix>
    void operator()(SizeInfo const& sizeInfo, Prefix prefix)
    {
      Impl::resizeVector(vector_, BlockingTag{}, sizeInfo, prefix);
    }

    template <class SizeInfo>
    void operator()(SizeInfo const& sizeInfo)
    {
      typename SizeInfo::SizePrefix prefix;
      prefix.resize(0);
      (*this)(sizeInfo, prefix);
    }

  private:
    Vector& vector_;
  };


  namespace Impl
  {
    template <class M, class RowBlocking, class ColBlocking,
              class RowSize, class RowPrefix, class ColSize, class ColPrefix>
    static void resizeMatrix(M& mat, RowBlocking br, ColBlocking bc,
                             RowSize const& rows, RowPrefix rowPrefix,
                             ColSize const& cols, ColPrefix colPrefix)
    {
      using MatrixResizable = Std::bool_constant<Concept::MatrixResizable<M>()>;
      using RowBlocked = Std::bool_constant<Concept::Blocked<RowBlocking>()>;
      using ColBlocked = Std::bool_constant<Concept::Blocked<ColBlocking>()>;

      auto num_rows = rows.size(rowPrefix);
      auto num_cols = cols.size(colPrefix);

      Hybrid::ifElse(MatrixResizable{},
        [&](auto id) { resize(id(mat), num_rows, num_cols); },
        [&](auto id) {
          if (RowBlocked::value && !ColBlocked::value) num_cols = 1;
          else if (!RowBlocked::value && ColBlocked::value) num_rows = 1;

          assert(hybridNumRows(mat) == num_rows && hybridNumCols(mat) == num_cols && "Wrong (static) dimension.");
        });

      Hybrid::ifElse(Std::bool_constant<RowBlocked::value || ColBlocked::value>{}, [&](auto id)
      {
        if (RowBlocked::value) { rowPrefix.push_back(0); }
        if (ColBlocked::value) { colPrefix.push_back(0); }
        Hybrid::forEach(range(hybridNumRows(mat)), [&](auto const _i) {
          if (RowBlocked::value) { rowPrefix.back() = _i; }
          Hybrid::forEach(range(hybridNumCols(mat)), [&](auto const _j) {
            if (ColBlocked::value) { colPrefix.back() = _j; }
            resizeMatrix(id(mat)[_i][_j], br[_i], bc[_j], rows, rowPrefix, cols, colPrefix);
          });
        });
      });
    }

  } // end namespace Impl


  template <class Matrix, class RowBlocking, class ColBlocking>
  class HierarchicMatrixResizer
  {
  public:
    HierarchicMatrixResizer(Matrix& matrix)
      : matrix_(matrix)
    {}

    template <class RowSize, class RowPrefix, class ColSize, class ColPrefix>
    void operator()(RowSize const& rows, RowPrefix rowPrefix,
                    ColSize const& cols, ColPrefix colPrefix)
    {
      Impl::resizeMatrix(matrix_, RowBlocking{}, ColBlocking{}, rows, rowPrefix, cols, colPrefix);
    }

    template <class RowSize, class ColSize>
    void operator()(RowSize const& rows, ColSize const& cols)
    {
      typename RowSize::SizePrefix rowPrefix;
      typename ColSize::SizePrefix colPrefix;
      rowPrefix.resize(0);
      colPrefix.resize(0);

      (*this)(rows, rowPrefix, cols, colPrefix);
    }

  private:
    Matrix& matrix_;
  };

} // end namespace Functions
} // end namespace Dune
