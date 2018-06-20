// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH

#include <set>
#include <algorithm>

#include <dune/common/test/testsuite.hh>
#include <dune/common/concept.hh>

#include <dune/geometry/quadraturerules.hh>

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



/*
 * Check if basis functions are continuous across faces.
 * Continuity is checked by evaluation at a set of quadrature points
 * from a quadrature rule of given order.
 * If two basis functions (on neighboring elements) share the same
 * global index, their values at the quadrature points (located on
 * their intersection) should coincide up to the given tolerance.
 *
 * If e basis function only appears on one side of the intersection,
 * it should be zero on the intersection.
 */
template<class Basis>
Dune::TestSuite checkBasisContinuity(const Basis& basis, std::size_t order = 5, double tol = 1e-10)
{
  Dune::TestSuite test("Global continuity check of basis functions");


  auto localView = basis.localView();
  auto neighborLocalView = basis.localView();

  auto norm = [](const auto& x) {
    return x.infinity_norm();
  };

  auto dist = [norm](const auto& x, const auto& y) {
    auto diff = x;
    diff -= y;
    return norm(diff);
  };


  for (const auto& e : elements(basis.gridView()))
  {
    localView.bind(e);
    for(const auto& intersection : intersections(basis.gridView(), e))
    {
      if (intersection.neighbor())
      {
        auto quadRule = Dune::QuadratureRules<double, Basis::GridView::dimension-1>::rule(intersection.type(), order);

        neighborLocalView.bind(intersection.outside());

        Dune::TypeTree::forEachLeafNode(localView.tree(), [&](const auto& node, auto&& treePath) {
          const auto& neighborNode = Dune::TypeTree::child(neighborLocalView.tree(), treePath);

          using Range = typename std::decay_t<decltype(node)>::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
          std::vector<std::vector<Range>> values;
          std::vector<std::vector<Range>> neighborValues;

          values.resize(quadRule.size());
          neighborValues.resize(quadRule.size());
          for(std::size_t k=0; k<quadRule.size(); ++k)
          {
            auto pointInElement = intersection.geometryInInside().global(quadRule[k].position());
            auto pointInNeighbor = intersection.geometryInOutside().global(quadRule[k].position());
            node.finiteElement().localBasis().evaluateFunction(pointInElement, values[k]);
            neighborNode.finiteElement().localBasis().evaluateFunction(pointInNeighbor, neighborValues[k]);
          }

          for(std::size_t i=0; i<node.size(); ++i)
          {
            bool foundInNeighbor = false;
            double maxJump = 0.0;
            for(std::size_t j=0; j<neighborNode.size(); ++j)
            {
              if (localView.index(node.localIndex(i)) == neighborLocalView.index(neighborNode.localIndex(j)))
              {
                // Basis function should only appear once in the neighbor element.
                test.check(foundInNeighbor==false)
                  << "Basis function " << localView.index(node.localIndex(i))
                  << " appears twice in element "
                  << neighborLocalView.element().type() << "#" << basis.gridView().indexSet().index(neighborLocalView.element()) ;

                // If basis function appears in neighbor element, then the
                // jump should be (numerically) zero across the intersection.
                for(std::size_t k=0; k<quadRule.size(); ++k)
                  maxJump = std::max(maxJump, dist(values[k][i], neighborValues[k][j]));
              }
            }
            // If basis function does not appear in neighbor element, then it
            // should be (numerically) zero on the intersection.
            if (not(foundInNeighbor))
              for(std::size_t k=0; k<quadRule.size(); ++k)
                maxJump = std::max(maxJump, norm(values[k][i]));
            test.check(maxJump < tol)
              << "Basis function " << localView.index(node.localIndex(i))
              << " is discontinuous across intersection of elements "
              << localView.element().type() << "#" << basis.gridView().indexSet().index(localView.element()) << " and "
              << neighborLocalView.element().type() << "#" << basis.gridView().indexSet().index(neighborLocalView.element());
              // << " since maximal jump " << maxJump << " exceeds tolerance";
          }
        });
      }
    }
  }
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
