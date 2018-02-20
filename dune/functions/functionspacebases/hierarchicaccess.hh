// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICACCESS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICACCESS_HH

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <dune/functions/common/access.hh>
#include <dune/functions/common/hybridsize.hh>
#include <dune/functions/common/indexaccess.hh>
#include <dune/functions/functionspacebases/blocking.hh>

namespace Dune { namespace Functions
{
  namespace Impl
  {
    // Provide operator[] index-access for containers
    template <class V, class BlockingTag, class Index, class F,
      std::enable_if_t<Concept::VectorAccessible<V,Index>(), int> = 0>
    decltype(auto) subVectorAccess(V&& vec, BlockingTag b, const Index& i, F&& f)
    {
      return f(access(vec,i), b[i]);
    }

    // Provide operator[] integral_constant-access for containers
    template <class V, class BlockingTag, class Index, class F,
      std::enable_if_t<not Concept::VectorAccessible<V,Index>(), int> = 0>
    decltype(auto) subVectorAccess(V&& vec, BlockingTag b, const Index& i, F&& f)
    {
      using Size = decltype(hybridSize(vec));
      return Hybrid::switchCases(std::make_index_sequence<Size::value>(), i,
        [&](const auto ii) -> decltype(auto) {
          return f(access(vec,ii), b[ii]);
        }, [&]() -> decltype(auto) {
          return f(access(vec,Indices::_0), b[0]);
        });
    }


    // Provide operator[][] index-access for containers
    template <class M, class BR, class BC, class I, class J, class F,
      std::enable_if_t<Concept::MatrixAccessible<M,I,J>() || Concept::isCallable<M,I,J>(), int> = 0>
    decltype(auto) subMatrixAccess(M&& mat, BR br, BC bc, I const& i, J const& j, F&& f)
    {
      return f(access(mat,i,j), br[i], bc[j]);
    }

    // Provide operator[][] integral_constant-access for containers
    template <class M, class BR, class BC, class I, class J, class F,
      std::enable_if_t<not Concept::MatrixAccessible<M,I,J>() && not Concept::isCallable<M,I,J>(), int> = 0>
    decltype(auto) subMatrixAccess(M&& mat, BR br, BC bc, I const& i, J const& j, F&& f)
    {
      using NumRows = decltype(hybridNumRows(mat));
      using NumCols = decltype(hybridNumCols(mat));
      return Hybrid::switchCases(std::make_index_sequence<NumRows::value>(), i,
        [&](const auto ii) -> decltype(auto) {
          return Hybrid::switchCases(std::make_index_sequence<NumCols::value>(), j,
            [&](const auto jj) -> decltype(auto) {
              return f(access(mat,ii,jj), br[ii], bc[jj]);
            }, [&]() -> decltype(auto) {
              // fall-back implementation
              return f(access(mat,ii,0), br[ii], bc[0]);
            });
        }, [&]() -> decltype(auto) {
          // fall-back implementation
          return f(access(mat,Indices::_0,Indices::_0), br[0], bc[0]);
        });
    }


    /// \brief Helper-class for the recursive operator[] calls to a vector.
    /**
     * Used as functor in \ref subVectorAccess and walks through the entries
     * of the multiindices using a \ref ShiftedMultiIndex helper.
     **/
    template <class MultiIndex>
    struct MultiIndexVectorResolver
    {
      MultiIndexVectorResolver(MultiIndex const& index)
        : index_(index)
      {}

      template <class C>
      decltype(auto) operator()(C&& c, Blocking::Flat) const
      {
        return c[index_[Indices::_0]];
      }

      template <class C, class BlockingTag>
      decltype(auto) operator()(C&& c, BlockingTag b) const
      {
        auto&& subIndex = shiftedMultiIndex(index_);
        auto&& subIndexResolver = MultiIndexVectorResolver<decltype(subIndex)>(subIndex);
        return subVectorAccess(c, b, index_[Indices::_0], subIndexResolver);
      }

      MultiIndex const& index_;
    };


    /// \brief Helper-class for the recursive operator[][] calls to a matrix.
    /**
     * Used as functor in \ref subMatrixAccess and walks through the entries
     * of the multiindices using a \ref ShiftedMultiIndex helper.
     **/
    template <class RowIndex, class ColIndex>
    struct MultiIndexMatrixResolver
    {
      MultiIndexMatrixResolver(RowIndex const& row, ColIndex const& col)
        : row_(row)
        , col_(col)
      {}

      template <class M>
      decltype(auto) operator()(M&& mat, Blocking::Flat, Blocking::Flat) const
      {
        return access(mat, row_[Indices::_0], col_[Indices::_0]);
      }

      template <class M, class RowBlocking>
      decltype(auto) operator()(M&& mat, RowBlocking br, Blocking::Flat bc) const
      {
        auto&& subRow = shiftedMultiIndex(row_);
        auto&& subIndexResolver = MultiIndexMatrixResolver<decltype(subRow),ColIndex>(subRow,col_);
        return subMatrixAccess(mat, br, bc, row_[Indices::_0], Indices::_0, subIndexResolver);
      }

      template <class M, class ColBlocking>
      decltype(auto) operator()(M&& mat, Blocking::Flat br, ColBlocking bc) const
      {
        auto&& subCol = shiftedMultiIndex(col_);
        auto&& subIndexResolver = MultiIndexMatrixResolver<RowIndex, decltype(subCol)>(row_,subCol);
        return subMatrixAccess(mat, br, bc, Indices::_0, col_[Indices::_0], subIndexResolver);
      }

      template <class M, class RowBlocking, class ColBlocking>
      decltype(auto) operator()(M&& mat, RowBlocking br, ColBlocking bc) const
      {
        auto&& subRow = shiftedMultiIndex(row_);
        auto&& subCol = shiftedMultiIndex(col_);
        auto&& subIndexResolver = MultiIndexMatrixResolver<decltype(subRow),decltype(subCol)>(subRow,subCol);
        return subMatrixAccess(mat, br, bc, row_[Indices::_0], col_[Indices::_0], subIndexResolver);
      }

      RowIndex const& row_;
      ColIndex const& col_;
    };

  } // end namespace Impl


  /// \brief Provide multiindex access to a vector for a given blocking structure
  /**
   * This provides access to a nested container by given
   * multiindex. Internally this is resolved by recusive
   * operator[]-calls with static or dynamic indices.
   * Because this recursion must be terminated using a
   * compile-time criterion, the blocking structure of
   * the container must explicitly be provided. The
   * recursion will terminate once the flat structure
   * \ref Blocking::Flat is reached.
   *
   * \param vec   The (hierarchic) container
   * \param b     The blocking structure of the container. See namespace \ref Blocking
   * \param index The multiindex of the entry you want to access.
   **/
  template <class Vector, class BlockingTag, class MultiIndex>
  decltype(auto) vectorAccess(Vector&& vec, BlockingTag const& b, const MultiIndex& index)
  {
    Impl::MultiIndexVectorResolver<MultiIndex> multiIndexResolver(index);
    return multiIndexResolver(vec, b);
  }

  /// \brief Provide multiindex access to a matrix for a given row/col blocking structure
  /**
   * This provides access to a nested matrix by given row/column
   * multiindices. Internally this is resolved by recusive
   * operator[][]-calls with static or dynamic indices.
   * Because this recursion must be terminated using a
   * compile-time criterion, the blocking structure of
   * the rows and columns must explicitly be provided. The
   * recursion will terminate once the flat structures
   * \ref Blocking::Flat are reached.
   *
   * \param vec   The (hierarchic) matrix
   * \param br    The blocking structure of the rows. See namespace \ref Blocking
   * \param bc    The blocking structure of the columns. See namespace \ref Blocking
   * \param row   The row-multiindex of the entry you want to access.
   * \param col   The column-multiindex of the entry you want to access.
   **/
  template <class Matrix, class RowBlocking, class ColBlocking, class RowIndex, class ColIndex>
  decltype(auto) matrixAccess(Matrix&& mat, RowBlocking const& br, ColBlocking const& bc,
                              RowIndex const& row, ColIndex const& col)
  {
    Impl::MultiIndexMatrixResolver<RowIndex, ColIndex> multiIndexResolver(row, col);
    return multiIndexResolver(mat, br, bc);
  }

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICACCESS_HH
