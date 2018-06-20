// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH

#include <set>
#include <algorithm>
#include <string>
#include <sstream>

#include <dune/common/test/testsuite.hh>
#include <dune/common/concept.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/functions/functionspacebases/concepts.hh>



/*
 * Get string identifier of element
 */
template<class Element, class GridView>
std::string elementStr(const Element& element, const GridView& gridView)
{
  std::stringstream s;
  s << element.type() << "#" << gridView.indexSet().index(element);
  return s.str();
}

/*
 * Check if two multi-indices are consecutive.
 * This is a used by checkBasisIndexTreeConsistency()
 */
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



/*
 * Check if given set of multi-indices is consistent, i.e.,
 * if it induces a consistent ordered tree. This is used
 * by checkBasisIndices()
 */
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



/*
 * Check consistency of basis.size(prefix)
 */
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



/*
 * Check indices of basis:
 * - First store the whole index tree in a set
 * - Check if this corresponds to a consistent index tree
 * - Check if index tree is consistent with basis.size(prefix) and basis.dimension()
 */
template<class Basis>
Dune::TestSuite checkBasisIndices(const Basis& basis)
{
  Dune::TestSuite test("basis index check");

  using MultiIndex = typename Basis::MultiIndex;

  static_assert(Dune::IsIndexable<MultiIndex>(), "MultiIndex must support operator[]");

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
      for(auto mi: multiIndex)
        test.check(mi>=0)
          << "Global multi-index containes negative entry for shape function " << i
          << " in element " << elementStr(localView.element(), basis.gridView());
      multiIndexSet.insert(multiIndex);
    }
  }

  test.subTest(checkBasisIndexTreeConsistency(multiIndexSet));
  test.subTest(checkBasisSizeConsistency(basis, multiIndexSet));
  test.check(basis.dimension() == multiIndexSet.size())
    << "basis.dimension() does not equal the total number of basis functions.";

  return test;
}



/*
 * Check if shape functions are not constant zero.
 * This is called by checkLocalView().
 */
template<class LocalFiniteElement>
Dune::TestSuite checkNonZeroShapeFunctions(const LocalFiniteElement& fe, std::size_t order = 5, double tol = 1e-10)
{
  Dune::TestSuite test;
  static const int dimension = LocalFiniteElement::Traits::LocalBasisType::Traits::dimDomain;

  auto quadRule = Dune::QuadratureRules<double, dimension>::rule(fe.type(), order);

  std::vector<typename LocalFiniteElement::Traits::LocalBasisType::Traits::RangeType> values;
  std::vector<bool> isNonZero;
  isNonZero.resize(fe.size(), false);
  for (const auto& qp : quadRule)
  {
    fe.localBasis().evaluateFunction(qp.position(), values);
    for(std::size_t i=0; i<fe.size(); ++i)
      isNonZero[i] = (isNonZero[i] or (values[i].infinity_norm() > tol));
  }
  for(std::size_t i=0; i<fe.size(); ++i)
    test.check(isNonZero[i])
      << "Found a constant zero basis function";
  return test;
}



/*
 * Check localView. This especially checks for
 * consistency of local indices and local size.
 */
template<class Basis, class LocalView>
Dune::TestSuite checkLocalView(const Basis& basis, const LocalView& localView)
{
  Dune::TestSuite test(std::string("LocalView on ") + elementStr(localView.element(), basis.gridView()));

  test.check(localView.size() <= localView.maxSize(), "localView.size() check")
    << "localView.size() is " << localView.size() << " but localView.maxSize() is " << localView.maxSize();

  // Count all local indices appearing in the tree.
  std::vector<std::size_t> localIndices;
  localIndices.resize(localView.size(), 0);
  Dune::TypeTree::forEachLeafNode(localView.tree(), [&](const auto& node, auto&& treePath) {
    test.check(node.size() == node.finiteElement().size())
      << "Size of leaf node and finite element are different.";
    for(std::size_t i=0; i<node.size(); ++i)
    {
      test.check(node.localIndex(i) < localView.size())
        << "Local index exceeds localView.size().";
      if (node.localIndex(i) < localView.size())
        ++(localIndices[node.localIndex(i)]);
    }
  });

  // Check if each local index appears exactly once.
  for(std::size_t i=0; i<localView.size(); ++i)
  {
    if (localIndices[i])
    test.check(localIndices[i]>=1)
      << "Local index " << i << " did not appear";
    test.check(localIndices[i]<=1)
      << "Local index " << i << " appears multiple times";
  }

  // Check if all basis functions are non-constant.
  Dune::TypeTree::forEachLeafNode(localView.tree(), [&](const auto& node, auto&& treePath) {
    test.subTest(checkNonConstantShapeFunctions(node.finiteElement()));
  });

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
 * If a basis function only appears on one side of the intersection,
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
                  << " appears twice in element " << elementStr(neighborLocalView.element(), basis.gridView());
                foundInNeighbor = true;
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
              << elementStr(localView.element(), basis.gridView())
              << " and " << elementStr(neighborLocalView.element(), basis.gridView());
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

  // Check if basis models the GlobalBasis concept.
  test.check(Dune::models<Dune::Functions::Concept::GlobalBasis<GridView>, Basis>(), "global basis concept check")
    << "type passed to checkBasis() does not model the GlobalBasis concept";

  // Perform all local tests.
  auto localView = basis.localView();
  for (const auto& e : elements(basis.gridView()))
  {
    localView.bind(e);
    test.subTest(checkLocalView(basis, localView));
  }

  // Perform global index tests.
  test.subTest(checkBasisIndices(basis));

  return test;
}

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
