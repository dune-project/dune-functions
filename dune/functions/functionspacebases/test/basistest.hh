// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH

#include <set>
#include <algorithm>
#include <string>
#include <sstream>
#include <map>
#include <optional>

#include <dune/common/test/testsuite.hh>
#include <dune/common/concept.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/common/typetree/childaccess.hh>
#include <dune/common/typetree/traversal.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/test/enabledifferentiabilitycheck.hh>
#include <dune/functions/functionspacebases/test/testboundlocalfe.hh>

struct CheckBasisFlag {};
struct AllowZeroBasisFunctions {};

// Enable checks that compare evaluateJacobian/partial methods up to order diffOrder with a
// finite difference approximation, check duality, and representation of constants on each
// element in the gridView
template<int i = 1>
struct CheckLocalFiniteElementFlag
{
  static constexpr int diffOrder = i;
};

template<class T, class... S>
struct IsContained : public std::disjunction<std::is_same<T,S>...>
{};

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

  // Based on the index tree, build a map that contains all possible
  // prefixes and maps each prefix to the size (of the subsequent digit).
  using Prefix = typename Basis::SizePrefix;
  auto prefixSet = std::map<Prefix, std::size_t>();
  for(const auto& index : multiIndexSet)
  {
    auto prefix  = Prefix();
    for (const auto& i: index)
    {
      prefixSet[prefix] = std::max(prefixSet[prefix], i+1);
      prefix.push_back(i);
    }
    prefixSet[prefix] = 0;
  }

  // Now check for all prefixes, if the size computed from the
  // index tree is consistent with basis.size(prefix).
  for(const auto& [prefix, size] : prefixSet)
  {
    auto prefixSize = basis.size(prefix);
    test.check(prefixSize == size, "basis.size(prefix) check")
      << "basis.size(" << prefix << ")=" << prefixSize << ", but should be " << size;
  }

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

    for (decltype(localView.size()) i=0; i< localView.size(); ++i)
    {
      auto multiIndex = localView.index(i);
      for(auto mi: multiIndex)
        test.check(mi>=0)
          << "Global multi-index contains negative entry for shape function " << i
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

  /**
   * Check that finite element returned by localView fulfills properties of local finite element
   * This test corresponds to a dune-localfunctions test, but for a bound, i.e. possibly transformed, FE
   *  This is called by checkLocalView()
   */

  template <int diffOrder, class LocalFiniteElement, class Element>
  Dune::TestSuite checkLocalFiniteElement(const LocalFiniteElement &fe, Element const &element)
  {
    Dune::TestSuite test("Local Finite Element Test");
    test.check(testBoundFE<diffOrder>(fe, element));
    return test;
  }
/*
 * Check localView. This especially checks for
 * consistency of local indices and local size.
 */
template<class Basis, class LocalView, class... Flags>
Dune::TestSuite checkLocalView(const Basis& basis, const LocalView& localView, Flags... flags)
{
  if (not localView.bound())
  {
    Dune::TestSuite test(std::string("Unbound LocalView"));
    Dune::TypeTree::forEachNode(localView.tree(), [&](const auto& node, auto&& treePath) {
      test.check(node.size()==0)
          << "Node of unbound LocalView has non-zero size " << node.size();
    });
    test.check(localView.size()==0)
        << "Unbound LocalView has non-zero size " << localView.size();
    return test;
  }

  Dune::TestSuite test(std::string("LocalView on ") + elementStr(localView.element(), basis.gridView()));

  test.check(localView.size() <= localView.maxSize(), "localView.size() check")
    << "localView.size() is " << localView.size() << " but localView.maxSize() is " << localView.maxSize();

  // Count all local indices appearing in the tree.
  std::vector<std::size_t> localIndices;
  localIndices.resize(localView.size(), 0);
  Dune::TypeTree::forEachLeafNode(localView.tree(), [&](const auto& node, auto&& treePath) {
    if (node.empty()) return;
    test.check(node.size() == node.finiteElement().size())
      << "Size of leaf node and finite element are different.";
    for(std::size_t i=0; i<node.size(); ++i)
    {
      test.check(node.localIndex(i) < localView.size())
        << "Local index exceeds localView.size() for leaf node " << treePath << ". "
        << node.localIndex(i) << " <? " <<  localView.size();
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

  // Check that all basis functions are non-zero.
  if (not IsContained<AllowZeroBasisFunctions, Flags...>::value)
  {
    Dune::TypeTree::forEachLeafNode(localView.tree(), [&](const auto& node, auto&& treePath) {
      if (node.empty()) return;
      test.subTest(checkNonZeroShapeFunctions(node.finiteElement()));
    });
  }
  if constexpr(IsContained<CheckLocalFiniteElementFlag<0>, Flags...>::value)
  {
    auto e = localView.element();
    Dune::TypeTree::forEachLeafNode(
        localView.tree(), [&](const auto &node, [[maybe_unused]] auto &&treePath)
        { test.subTest(checkLocalFiniteElement<0>(node.finiteElement(), e)); });
  } else if constexpr (IsContained<CheckLocalFiniteElementFlag<1>,
                              Flags...>::value) {
      auto e = localView.element();
      Dune::TypeTree::forEachLeafNode(
          localView.tree(), [&](const auto &node, [[maybe_unused]] auto &&treePath)
          { test.subTest(checkLocalFiniteElement<1>(node.finiteElement(), e)); });
  } else if constexpr (IsContained<CheckLocalFiniteElementFlag<2>,
                                    Flags...>::value) {
      auto e = localView.element();
      Dune::TypeTree::forEachLeafNode(
          localView.tree(), [&](const auto &node, [[maybe_unused]] auto &&treePath)
          { test.subTest(checkLocalFiniteElement<2>(node.finiteElement(), e)); });
  }

  // Check if copies of the local view are independent.
  // To this end we create a copy, bind it to another element
  // and check if the original LocalView still behaves correct.
  // The check is done by first creating a fresh LocalView bound
  // to the same element as reference.
  {
    // Create a new local view as reference
    auto localView_reference = basis.localView();
    localView_reference.bind(localView.element());

    // Create and modify a copy of the local view
    auto localView_copy = localView;
    localView_copy.bind(*basis.gridView().template begin<0>());

    Dune::TypeTree::forEachNode(localView.tree(), [&](const auto& node, auto&& treePath) {
      const auto& node_reference = Dune::TypeTree::child(localView_reference.tree(), treePath);
      const auto& node_copy = Dune::TypeTree::child(localView_copy.tree(), treePath);
      test.check(&node != &node_copy)
        << "LocalView and its copy share node " << treePath;
      test.check(node.size() == node_reference.size())
        << "Size of node changed by modifying a copy of the LocalView.";
      if (node.empty()) return;
      for(std::size_t i=0; i<node.size(); ++i)
      {
        test.check(node.localIndex(i) == node_reference.localIndex(i))
          << "Local index of node changed by modifying a copy of the LocalView.";
        test.check(localView.index(node.localIndex(i)) == localView_reference.index(node_reference.localIndex(i)))
          << "Global index of node changed by modifying a copy of the LocalView.";
      }
    });
  }

  return test;
}


// Flag to enable a local continuity check for checking strong
// continuity across an intersection within checkBasisContinuity().
//
// For each inside basis function this will compute the jump against
// zero or the corresponding inside basis function. The latter is then
// checked for being (up to a tolerance) zero on a set of quadrature points.
struct EnableContinuityCheck
{
  std::size_t order_ = 5;
  double tol_ = 1e-10;

  template<class JumpEvaluator>
  auto localJumpContinuityCheck(const JumpEvaluator& jumpEvaluator, std::size_t order, double tol) const
  {
    return [=](const auto& intersection, const auto& treePath, const auto& insideNode, const auto& outsideNode, const auto& insideToOutside) {
      using Intersection = std::decay_t<decltype(intersection)>;
      using Node = std::decay_t<decltype(insideNode)>;

      std::vector<int> isContinuous(insideNode.size(), true);
      const auto& quadRule = Dune::QuadratureRules<double, Intersection::mydimension>::rule(intersection.type(), order);

      using Range = typename Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
      std::vector<std::vector<Range>> values;
      std::vector<std::vector<Range>> neighborValues;

      // Evaluate inside and outside basis functions.
      values.resize(quadRule.size());
      neighborValues.resize(quadRule.size());
      for(std::size_t k=0; k<quadRule.size(); ++k)
      {
        auto pointInElement = intersection.geometryInInside().global(quadRule[k].position());
        auto pointInNeighbor = intersection.geometryInOutside().global(quadRule[k].position());
        insideNode.finiteElement().localBasis().evaluateFunction(pointInElement, values[k]);
        outsideNode.finiteElement().localBasis().evaluateFunction(pointInNeighbor, neighborValues[k]);
      }

      // Check jump against outside basis function or zero.
      for(std::size_t i=0; i<insideNode.size(); ++i)
      {
        for(std::size_t k=0; k<quadRule.size(); ++k)
        {
          auto jump = values[k][i];
          if (insideToOutside[i].has_value())
            jump -= neighborValues[k][insideToOutside[i].value()];
          isContinuous[i] = isContinuous[i] and (jumpEvaluator(jump, intersection, quadRule[k].position()) < tol);
        }
      }
      return isContinuous;
    };
  }

  auto localContinuityCheck() const {
    auto jumpNorm = [](auto&&jump, auto&& intersection, auto&& x) -> double {
      return jump.infinity_norm();
    };
    return localJumpContinuityCheck(jumpNorm, order_, tol_);
  }
};

// Flag to enable a local normal-continuity check for checking strong
// continuity across an intersection within checkBasisContinuity().
//
// For each inside basis function this will compute the normal jump against
// zero or the corresponding inside basis function. The latter is then
// checked for being (up to a tolerance) zero on a set of quadrature points.
struct EnableNormalContinuityCheck : public EnableContinuityCheck
{
  auto localContinuityCheck() const {
    auto normalJump = [](auto&&jump, auto&& intersection, auto&& x) -> double {
      return jump * intersection.unitOuterNormal(x);
    };
    return localJumpContinuityCheck(normalJump, order_, tol_);
  }
};

// Flag to enable a local tangential-continuity check for checking continuity
// of tangential parts of a vector-valued basis across an intersection
// within checkBasisContinuity().
//
// For each inside basis function this will compute the tangential jump against
// zero or the corresponding outside basis function. The jump is then
// checked for being (up to a tolerance) zero on a set of quadrature points.
struct EnableTangentialContinuityCheck : public EnableContinuityCheck
{
  auto localContinuityCheck() const {
    auto tangentialJumpNorm = [](auto&&jump, auto&& intersection, auto&& x) -> double {
      auto tangentialJump = jump - (jump * intersection.unitOuterNormal(x)) * intersection.unitOuterNormal(x);
      return tangentialJump.two_norm();
    };
    return localJumpContinuityCheck(tangentialJumpNorm, order_, tol_);
  }
};

// Flag to enable a center continuity check for checking continuity in the
// center of an intersection within checkBasisContinuity().
//
// For each inside basis function this will compute the jump against
// zero or the corresponding inside basis function. The latter is then
// checked for being (up to a tolerance) zero in the center of mass
// of the intersection.
struct EnableCenterContinuityCheck : public EnableContinuityCheck
{
  template<class JumpEvaluator>
  auto localJumpCenterContinuityCheck(const JumpEvaluator& jumpEvaluator, double tol) const
  {
    return [=](const auto& intersection, const auto& treePath, const auto& insideNode, const auto& outsideNode, const auto& insideToOutside) {
      using Node = std::decay_t<decltype(insideNode)>;
      using Range = typename Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;

      std::vector<int> isContinuous(insideNode.size(), true);
      std::vector<Range> insideValues;
      std::vector<Range> outsideValues;

      insideNode.finiteElement().localBasis().evaluateFunction(intersection.geometryInInside().center(), insideValues);
      outsideNode.finiteElement().localBasis().evaluateFunction(intersection.geometryInOutside().center(), outsideValues);

      auto centerLocal = intersection.geometry().local(intersection.geometry().center());

      // Check jump against outside basis function or zero.
      for(std::size_t i=0; i<insideNode.size(); ++i)
      {
          auto jump = insideValues[i];
          if (insideToOutside[i].has_value())
            jump -= outsideValues[insideToOutside[i].value()];
          isContinuous[i] = isContinuous[i] and (jumpEvaluator(jump, intersection, centerLocal) < tol);
      }
      return isContinuous;
    };
  }

  auto localContinuityCheck() const {
    auto jumpNorm = [](auto&&jump, auto&& intersection, auto&& x) -> double {
      return jump.infinity_norm();
    };
    return localJumpCenterContinuityCheck(jumpNorm, tol_);
  }
};

  // Flag to enable a vertex continuity check for checking continuity at the vertices
  // of an intersection within checkBasisContinuity().
  //
  // For each inside basis function this will compute the jump against
  // zero or the corresponding inside basis function. The latter is then
  // checked for being (up to a tolerance) zero in the vertices of mass
  // of the intersection.

  struct EnableVertexContinuityCheck: public EnableContinuityCheck
  {
    template <class JumpEvaluator>
    auto localJumpVertexContinuityCheck(const JumpEvaluator &jumpEvaluator, double tol) const
    {
      return [=](const auto &intersection, const auto &treePath, const auto &insideNode,
                 const auto &outsideNode, const auto &insideToOutside)
      {
        using Node = std::decay_t<decltype(insideNode)>;
        using Range = typename Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;

        std::vector<int> isContinuous(insideNode.size(), true);
        std::vector<Range> insideValues;
        std::vector<Range> outsideValues;

        std::vector<typename std::decay_t<decltype(intersection)>::LocalCoordinate> vertices;

        for (int i = 0; i < intersection.geometry().corners(); ++i)
          vertices.push_back(intersection.geometry().local(intersection.geometry().corner(i)));

        for (auto const &vertex : vertices)
        {

          insideNode.finiteElement().localBasis().evaluateFunction(
              intersection.geometryInInside().global(vertex), insideValues);
          outsideNode.finiteElement().localBasis().evaluateFunction(
              intersection.geometryInOutside().global(vertex), outsideValues);
          // Check jump against outside basis function or zero.
          for (std::size_t i = 0; i < insideNode.size(); ++i)
          {
            auto jump = insideValues[i];
            if (insideToOutside[i].has_value())
              jump -= outsideValues[insideToOutside[i].value()];
            isContinuous[i] = isContinuous[i] and (jumpEvaluator(jump, intersection, vertex) < tol);
          }
        }
        return isContinuous;
      };
    }

    auto localContinuityCheck() const
    {
      auto jumpNorm = [](auto &&jump, auto &&intersection, auto &&x) -> double
      { return jump.infinity_norm(); };
      return localJumpVertexContinuityCheck(jumpNorm, tol_);
    }
  };

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
template<class Basis, class LocalCheck>
Dune::TestSuite checkBasisContinuity(const Basis& basis, const LocalCheck& localCheck)
{
  Dune::TestSuite test("Global continuity check of basis functions");


  auto localView = basis.localView();
  auto neighborLocalView = basis.localView();

  for (const auto& e : elements(basis.gridView()))
  {
    localView.bind(e);
    for(const auto& intersection : intersections(basis.gridView(), e))
    {
      if (intersection.neighbor())
      {
        neighborLocalView.bind(intersection.outside());

        Dune::TypeTree::forEachLeafNode(localView.tree(), [&](const auto& insideNode, auto&& treePath) {
          if (insideNode.empty()) return;

          const auto& outsideNode = Dune::TypeTree::child(neighborLocalView.tree(), treePath);

          std::vector<std::optional<int>> insideToOutside;
          insideToOutside.resize(insideNode.size());

          // Map all inside DOFs to outside DOFs if possible
          for(std::size_t i=0; i<insideNode.size(); ++i)
          {
            for(std::size_t j=0; j<outsideNode.size(); ++j)
            {
              if (localView.index(insideNode.localIndex(i)) == neighborLocalView.index(outsideNode.localIndex(j)))
              {
                // Basis function should only appear once in the neighbor element.
                test.check(not insideToOutside[i].has_value())
                  << "Basis function " << localView.index(insideNode.localIndex(i))
                  << " appears twice in element " << elementStr(neighborLocalView.element(), basis.gridView());
                insideToOutside[i] = j;
              }
            }
          }

          // Apply continuity check on given intersection with given inside/outside DOF node pair.
          auto isContinuous = localCheck(intersection, treePath, insideNode, outsideNode, insideToOutside);

          for(std::size_t i=0; i<insideNode.size(); ++i)
          {
            test.check(isContinuous[i])
              << "Basis function " << localView.index(insideNode.localIndex(i))
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

template<class Basis, class... Flags>
Dune::TestSuite checkConstBasis(const Basis& basis, Flags... flags)
{
  Dune::TestSuite test("const basis check");

  using GridView = typename Basis::GridView;

  // Check if basis models the GlobalBasis concept.
  test.check(Dune::models<Dune::Functions::Concept::GlobalBasis<GridView>, Basis>(), "global basis concept check")
    << "type passed to checkBasis() does not model the GlobalBasis concept";

  // Perform all local tests.
  auto localView = basis.localView();
  test.check(not localView.bound())
    << "Unbound LocalView returns localView.bound()==true";
  test.subTest(checkLocalView(basis, localView, flags...));
  for (const auto& e : elements(basis.gridView()))
  {
    localView.bind(e);
    test.check(localView.bound())
      << "Bound LocalView returns localView.bound()==false";
    test.subTest(checkLocalView(basis, localView, flags...));
  }

  // Perform global index tests.
  test.subTest(checkBasisIndices(basis));

  // Perform continuity check.
  // First capture flags in a tuple in order to iterate.
  auto flagTuple = std::tie(flags...);
  Dune::Hybrid::forEach(flagTuple, [&](auto&& flag) {
    using Flag = std::decay_t<decltype(flag)>;
    if constexpr (std::is_base_of_v<EnableContinuityCheck, Flag>)
      test.subTest(checkBasisContinuity(basis, flag.localContinuityCheck()));
    else if constexpr (std::is_base_of_v<EnableDifferentiabilityCheck, Flag>)
      test.subTest(checkBasisDifferentiability(basis, flag));
  });

  return test;
}

template<class Basis, class... Flags>
Dune::TestSuite checkBasis(Basis& basis, Flags... flags)
{
  Dune::TestSuite test("basis check");

  // Perform tests for a constant basis
  test.subTest(checkConstBasis(basis,flags...));

  // Check copy-construction / copy-assignable of a basis
  {
    Basis copy(basis);
    test.subTest(checkConstBasis(copy,flags...));
    copy = basis;
    test.subTest(checkConstBasis(copy,flags...));
  }

  // Check update of gridView
  auto gridView = basis.gridView();
  basis.update(gridView);

  return test;
}




#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BASISTEST_HH
