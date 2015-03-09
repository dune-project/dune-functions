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



class InvalidRange
{};


template<class Signature>
struct DefaultDerivativeTraits
{
  typedef InvalidRange Range;
};

template<>
struct DefaultDerivativeTraits< double(double) >
{
  typedef double Range;
};

template<typename K, int n>
struct DefaultDerivativeTraits<K(FieldVector<K,n>)>
{
  typedef FieldVector<K,n> Range;
};

template<typename K, int n, int m>
struct DefaultDerivativeTraits<FieldVector<K,m>(FieldVector<K,n>)>
{
  typedef FieldMatrix<K,m,n> Range;
};

template<typename K, int n, int m>
struct DefaultDerivativeTraits<FieldMatrix<K,1,m>(FieldVector<K,n>)>
{
  typedef FieldMatrix<K,m,n> Range;
};


}} // namespace Dune::Functions


#endif // DUNE_FUNCTIONS_COMMON_DEFAULT_DERIVATIVE_TRAITS_HH
