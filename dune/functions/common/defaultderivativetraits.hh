// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DEFAULT_DERIVATIVE_TRAITS_HH
#define DUNE_FUNCTIONS_COMMON_DEFAULT_DERIVATIVE_TRAITS_HH

#include <type_traits>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune {
namespace Functions {



/**
 * \brief Dummy range class to be used if no proper type is available
 *
 * \ingroup FunctionUtility
 */
class InvalidRange
{};


/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * This class provides sensible defaults for the range
 * of derivatives of functions with some common \p Domain
 * and \p Range types.
 */
template<class Signature>
struct DefaultDerivativeTraits
{
  //! Range of derivative for function with given signature
  typedef InvalidRange Range;
};


/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * Specialization for Signature = double(double)
 */
template<>
struct DefaultDerivativeTraits< double(double) >
{
  //! \copydoc DefaultDerivativeTraits::Range
  typedef double Range;
};

/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * \tparam K Scalar range type
 *
 * Specialization for Signature = K(FieldVector<K,n>)
 */
template<typename K, int n>
struct DefaultDerivativeTraits<K(FieldVector<K,n>)>
{
  //! \copydoc DefaultDerivativeTraits::Range
  typedef FieldVector<K,n> Range;
};

/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * \tparam K Scalar range type
 *
 * Specialization for Signature = FieldVector<K,m>(FieldVector<K,n>)
 */
template<typename K, int n, int m>
struct DefaultDerivativeTraits<FieldVector<K,m>(FieldVector<K,n>)>
{
  //! \copydoc DefaultDerivativeTraits::Range
  typedef FieldMatrix<K,m,n> Range;
};

/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * \tparam K Scalar range type
 *
 * Specialization for Signature = FieldMatrix<K,1,m>(FieldVector<K,n>)
 */
template<typename K, int n, int m>
struct DefaultDerivativeTraits<FieldMatrix<K,1,m>(FieldVector<K,n>)>
{
  //! \copydoc DefaultDerivativeTraits::Range
  typedef FieldMatrix<K,m,n> Range;
};


}} // namespace Dune::Functions


#endif // DUNE_FUNCTIONS_COMMON_DEFAULT_DERIVATIVE_TRAITS_HH
