// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH

#include <set>
#include <algorithm>

#include <dune/common/test/testsuite.hh>
#include <dune/common/concept.hh>

#include <dune/functions/functionspacebases/concepts.hh>



// check if two multi-indices are consecutive
template<class MultiIndex>
bool multiIndicesConsecutive(const MultiIndex& a, const MultiIndex& b)
{
  std::size_t i = 0;

  // find largest common prefix
  for (; (i<a.size()) and (i<b.size()) and (a[i] == b[i]); ++i)
  {};

  // if b is exhausted but a is not, then b is a strict prefix of a and does not succeed a
  if ((i<a.size()) and (i==b.size()))
    return false;

  // if a and b are not exhausted, then the first non-common index must be an increment
  if ((i<a.size()) and (i<b.size()))
  {
    if (b[i] != a[i]+1)
      return false;
    ++i;
  }

  // if b is not exhausted, then the following indices should be zero
  if (i<b.size())
  {
    for (; i<b.size(); ++i)
    {
      if (b[i] != 0)
        return false;
    }
  }
  return true;
}



template<class MultiIndexSet>
Dune::TestSuite checkBasisIndexTreeConsistency(const MultiIndexSet& multiIndexSet)
{
  Dune::TestSuite test("index tree consistency check");

  using namespace Dune;

  auto it = multiIndexSet.begin();
  auto end = multiIndexSet.end();

  // get first multi-index
  auto lastMultiIndex = *it;

  // assert that index is non-empty
  test.require(lastMultiIndex.size()>0, "multi-index size check")
    << "empty multi-index found";

  // check if first multi-index is [0,...,0]
  for (decltype(lastMultiIndex.size()) i = 0; i<lastMultiIndex.size(); ++i)
  {
    test.require(lastMultiIndex[i] == 0, "smallest index check")
      << "smallest index contains non-zero entry " << lastMultiIndex[i] << " in position " << i;
  }

  ++it;
  for(; it != end; ++it)
  {
    auto multiIndex = *it;

    // assert that index is non-empty
    test.require(multiIndex.size()>0, "multi-index size check")
      << "empty multi-index found";

    // assert that indices are consecutive
    test.check(multiIndicesConsecutive(lastMultiIndex, multiIndex), "consecutive index check")
      << "multi-indices " << lastMultiIndex << " and " << multiIndex << " are subsequent but not consecutive";

    lastMultiIndex = multiIndex;
  }

  return test;
}



template<class Basis, class MultiIndexSet>
Dune::TestSuite checkBasisSizeConsistency(const Basis& basis, const MultiIndexSet& multiIndexSet)
{
  Dune::TestSuite test("index size consistency check");

  auto prefix = typename Basis::SizePrefix{};

  for(const auto& index : multiIndexSet)
  {
    prefix.clear();
    for (const auto& i: index)
    {
      // All indices i collected so far from the multi-index
      // refer to a non-empty multi-index subtree. Hence the
      // size must be nonzero and in fact strictly larger than
      // the next index.
      auto prefixSize = basis.size(prefix);
      test.require(prefixSize > i, "basis.size(prefix) subtree check")
        << "basis.size(" << prefix << ")=" << prefixSize << " but index " << index << " exists";

      // append next index from multi-index
      prefix.push_back(i);
    }
    auto prefixSize = basis.size(prefix);
    test.require(prefixSize == 0, "basis.size(prefix) leaf check")
      << "basis.size(" << prefix << ")=" << prefixSize << " but the prefix exists as index";
  }

  // ToDo: Add check that for basis.size(prefix)==n with i>0
  // there exist multi-indices of the form (prefix,0,...)...(prefix,n-1,...)

  return test;
}



template<class Basis>
Dune::TestSuite checkBasisIndices(const Basis& basis)
{
  Dune::TestSuite test("basis index check");

  using MultiIndex = typename Basis::MultiIndex;
  auto compare = [](const auto& a, const auto& b) {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
  };

  auto multiIndexSet = std::set<MultiIndex, decltype(compare)>{compare};

  auto localView = basis.localView();
  for (const auto& e : elements(basis.gridView()))
  {
    localView.bind(e);

    test.require(localView.size() <= localView.maxSize(), "localView.size() check")
      << "localView.size() is " << localView.size() << " but localView.maxSize() is " << localView.maxSize();

    for (decltype(localView.size()) i=0; i< localView.size(); ++i)
    {
      auto multiIndex = localView.index(i);
      multiIndexSet.insert(multiIndex);
    }
  }

  test.subTest(checkBasisIndexTreeConsistency(multiIndexSet));
  test.subTest(checkBasisSizeConsistency(basis, multiIndexSet));

  return test;
}



template<class Basis>
Dune::TestSuite checkBasis(const Basis& basis)
{
  Dune::TestSuite test("basis check");


  using GridView = typename Basis::GridView;

  test.check(Dune::models<Dune::Functions::Concept::GlobalBasis<GridView>, Basis>(), "global basis concept check")
    << "type passed to checkBasis() does not model the GlobalBasis concept";

  test.subTest(checkBasisIndices(basis));

  return test;
}

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
