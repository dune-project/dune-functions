// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH

#include <set>
#include <algorithm>

#include <dune/common/test/testsuite.hh>


template<class MultiIndexSet>
Dune::TestSuite checkIndexTreeConsistency(const MultiIndexSet& multiIndexSet)
{
  Dune::TestSuite test("index tree consistency check");

  auto it = multiIndexSet.begin();
  auto end = multiIndexSet.end();

  // check if two multi-indices are consecutive
  auto is_consecutive = [](auto&& a, auto&& b) {
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
  };

  // get first multi-index
  auto lastMultiIndex = *it;

  // assert that index is non-empty
  test.require(lastMultiIndex.size()>0, "multi-index size check")
    << "empty multi-index found";

  // check if first multi-index is [0,...,0]
  for (int i = 0; i<lastMultiIndex.size(); ++i)
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
    test.require(is_consecutive(lastMultiIndex, multiIndex), "consecutive index check")
      << "multi-indices " << lastMultiIndex << " and " << multiIndex << " are subsequent but not consecutive";
    lastMultiIndex = multiIndex;
  }

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
  auto localIndexSet = basis.localIndexSet();
  for (const auto& e : elements(basis.gridView()))
  {
    localView.bind(e);
    localIndexSet.bind(localView);

    test.require(localIndexSet.size() <= localView.maxSize(), "localIndexSet.size() check")
      << "localIndexSet.size() is " << localIndexSet.size() << " but localView.maxSize() is " << localView.maxSize();

    for(int i=0; i< localIndexSet.size(); ++i)
    {
      auto multiIndex = localIndexSet.index(i);
      multiIndexSet.insert(multiIndex);
    }
  }

  test.subTest(checkIndexTreeConsistency(multiIndexSet));

  return test;
}



template<class Basis>
Dune::TestSuite checkBasis(const Basis& basis)
{
  Dune::TestSuite test("basis check");

  test.subTest(checkBasisIndices(basis));

  return test;
}

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
