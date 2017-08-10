// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FUNCTIONRANGE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FUNCTIONRANGE_HH

#include <dune/common/fvector.hh>
#include <dune/common/concept.hh>
#include <dune/functions/functionspacebases/concepts.hh>

namespace Dune {
namespace Functions {

// Function range for vector valued functions
template<typename D, typename R, typename DHasIndexAccess = void, typename RHasIndexAccess = void>
struct FunctionRange {};



template<typename D, typename R>
struct FunctionRange<D, R,
                     typename std::enable_if< models<Concept::HasIndexAccess, D, int>() >::type,
                     typename std::enable_if< models<Concept::HasIndexAccess, R, int>() >::type>
{
  using Domain = D;
  using Range = R;
  using RangeField = typename Range::value_type;
  using JacobianRange = FieldVector< FieldVector<RangeField, Domain::dimension>, Range::dimension>;
};

template<typename D, typename R>
struct FunctionRange<D, R,
                     typename std::enable_if<     models<Concept::HasIndexAccess, D, int>() >::type,
                     typename std::enable_if< not models<Concept::HasIndexAccess, R, int>() >::type>
{
  using Domain = D;
  using Range = R;
  using RangeField = Range;
  using JacobianRange = FieldVector<RangeField, Domain::dimension>;
};

template<typename D, typename R>
struct FunctionRange<D, R,
                     typename std::enable_if< not models<Concept::HasIndexAccess, D, int>() >::type,
                     typename std::enable_if<     models<Concept::HasIndexAccess, R, int>() >::type>
{
  using Domain = D;
  using Range = R;
  using RangeField = typename Range::value_type;
  using JacobianRange = FieldVector<RangeField, Range::dimension>;
};

template<typename D, typename R>
struct FunctionRange<D, R,
                     typename std::enable_if< not models<Concept::HasIndexAccess, D, int>() >::type,
                     typename std::enable_if< not models<Concept::HasIndexAccess, R, int>() >::type>
{
  using Domain = D;
  using Range = R;
  using RangeField = Range;
  using JacobianRange = RangeField;
};


} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONRANGE_HH
