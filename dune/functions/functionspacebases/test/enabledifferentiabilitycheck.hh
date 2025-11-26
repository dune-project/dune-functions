#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_ENABLEDIFFERENTIABILITYCHECK_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_ENABLEDIFFERENTIABILITYCHECK_HH

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#include <optional>

#include <dune/common/test/testsuite.hh>
#include <dune/common/transpose.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/common/typetree/traversal.hh>

template <class Element, class GridView>
std::string elementStr(const Element &element, const GridView &gridView);

// Flag to enable a local continuity check for checking partial
// differentiability across an intersection within checkBasisDifferentiability().
//
// For each inside basis function this will compute the jump against
// zero or the corresponding inside basis function. The latter is then
// checked for being (up to a tolerance) zero on a set of quadrature points.
struct EnableDifferentiabilityCheck
{
  std::size_t order_ = 5;
  double tol_ = 1e-10;
  std::string checkLocation = "across";

  template <class JumpEvaluator, class QuadRuleFactoryMethod>
  auto localJumpDifferentiabilityCheck(const JumpEvaluator &jumpEvaluator,
                                        const QuadRuleFactoryMethod &quadRuleFactoryMethod,
                                        std::size_t order, double tol) const
  {
    return [=](const auto &intersection, const auto &treePath, const auto &insideNode,
                const auto &outsideNode, const auto &insideToOutside)
    {
      // using Intersection = std::decay_t<decltype(intersection)>;
      using Node = std::decay_t<decltype(insideNode)>;

      std::vector<int> isDifferentiable(insideNode.size(), true);
      const auto &quadRule = quadRuleFactoryMethod(intersection, order);

      // using Range = typename Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
      using JacobiRange =
          typename Node::FiniteElement::Traits::LocalBasisType::Traits::JacobianType;
      std::vector<std::vector<JacobiRange>> insideValues;
      std::vector<std::vector<JacobiRange>> outsideValues;

      // Evaluate inside and outside basis functions.
      insideValues.resize(quadRule.size());
      outsideValues.resize(quadRule.size());

      std::size_t insideNodeSize = insideNode.finiteElement().localBasis().size();
      std::size_t outsideNodeSize = outsideNode.finiteElement().localBasis().size();

      for (std::size_t k = 0; k < quadRule.size(); ++k)
      {
        std::vector<JacobiRange> insideLocalJacobians, outsideLocalJacobians;

        insideLocalJacobians.resize(insideNodeSize);
        outsideLocalJacobians.resize(outsideNodeSize);
        insideValues[k].resize(insideNodeSize);
        outsideValues[k].resize(outsideNodeSize);

        auto insidePoint = intersection.geometryInInside().global(quadRule[k].position());
        auto outsidePoint = intersection.geometryInOutside().global(quadRule[k].position());
        insideNode.finiteElement().localBasis().evaluateJacobian(insidePoint,
                                                                  insideLocalJacobians);
        outsideNode.finiteElement().localBasis().evaluateJacobian(outsidePoint,
                                                                  outsideLocalJacobians);
        auto insideJacobiInverseTransposed
            = intersection.inside().geometry().jacobianInverseTransposed(insidePoint);
        auto outsideJacobiInverseTransposed
            = intersection.outside().geometry().jacobianInverseTransposed(outsidePoint);
        for (std::size_t i = 0; i < insideNodeSize; ++i)
        {
          // TODO make this viable for localJac of type FieldVector (at least)
          insideValues[k][i] = insideLocalJacobians[i] * transpose(insideJacobiInverseTransposed);
          if (i < outsideNodeSize)
            outsideValues[k][i]
                = outsideLocalJacobians[i] * transpose(outsideJacobiInverseTransposed);
        }
      }

      // Check jump against outside basis function or zero.
      for (std::size_t i = 0; i < insideNode.size(); ++i)
      {
        for (std::size_t k = 0; k < quadRule.size(); ++k)
        {
          auto jump = insideValues[k][i];
          if (insideToOutside[i].has_value())
            jump -= outsideValues[k][insideToOutside[i].value()];
          isDifferentiable[i] = isDifferentiable[i]
                            and (jumpEvaluator(jump, intersection, quadRule[k].position()) < tol);
        }
      }
      for (std::size_t k = 0; k < quadRule.size(); ++k)
      {
        insideValues[k].clear();
        outsideValues[k].clear();
      }
      return isDifferentiable;
    };
  }

  auto localDifferentiabilityCheck() const
  {
    auto jumpNorm = [](auto &&jump, auto &&intersection, auto &&x) -> double
    { return jump.infinity_norm(); };

    auto quadRuleProvider = [](auto &&intersection, auto &&order)
    {
      return Dune::QuadratureRules<double, std::decay_t<decltype(intersection)>::mydimension>::rule(intersection.type(),
                                                                            order);
    };

    return localJumpDifferentiabilityCheck(jumpNorm, quadRuleProvider, order_, tol_);
  }
};

struct EnableVertexDifferentiabilityCheck: public EnableDifferentiabilityCheck
{
  std::size_t order_ = 5;
  double tol_ = 1e-10;
  std::string checkLocation = "at vertices of";

  auto localDifferentiabilityCheck() const
  {
    auto jumpNorm = [](auto &&jump, auto &&intersection, auto &&x) -> double
    { return jump.infinity_norm(); };

    auto quadRuleProvider = [](auto &&intersection, auto &&order)
    {
      using Pt = Dune::QuadraturePoint<double, std::decay_t<decltype(intersection)>::mydimension>;
      std::vector<Pt> quadRule;
      for (int i = 0; i < intersection.geometry().corners(); ++i)
        quadRule.push_back(
            Pt{intersection.geometry().local(intersection.geometry().corner(i)), 1.});
      return quadRule;
    };

    return localJumpDifferentiabilityCheck(jumpNorm, quadRuleProvider, order_, tol_);
  }
};

struct EnableNormalDifferentiabilityAtMidpointsCheck: public EnableDifferentiabilityCheck
{
  std::size_t order_ = 5;
  double tol_ = 1e-10;
  std::string checkLocation = "at edge midpoint of";

  auto localDifferentiabilityCheck() const
  {
    auto jumpNorm = [](auto &&jump, auto &&intersection, auto &&x) -> double
    {
      Dune::FieldVector<double, 1> res;
      jump.mv(intersection.unitOuterNormal(x), res);
      return res.infinity_norm();
    };

    auto quadRuleProvider = [](auto &&intersection, auto &&order)
    {
      using Pt = Dune::QuadraturePoint<double, std::decay_t<decltype(intersection)>::mydimension>;
      std::vector<Pt> quadRule;

      quadRule.push_back(Pt{intersection.geometry().local(intersection.geometry().center()), 1.});
      return quadRule;
    };

    return localJumpDifferentiabilityCheck(jumpNorm, quadRuleProvider, order_, tol_);
  }
};

/*
  * Check if basis functions are differentiable across faces.
  * Differeniability is checked by evaluation at a set of quadrature points
  * from a quadrature rule of given order, or at a given set of points, like vertices/ center, etc.
  * If two basis functions (on neighboring elements) share the same
  * global index, their derivatives at the quadrature points (located on
  * their intersection) should coincide up to the given tolerance.
  *
  * If a basis function only appears on one side of the intersection,
  * it should be zero on the intersection.
  */
template <class Basis, class Flag>
Dune::TestSuite checkBasisDifferentiability(const Basis &basis, const Flag &flag)
{
  Dune::TestSuite test("Global differentiability check of basis functions");

  auto const &localCheck = flag.localDifferentiabilityCheck();

  auto localView = basis.localView();
  auto neighborLocalView = basis.localView();

  for (const auto &e : elements(basis.gridView()))
  {
    localView.bind(e);
    for (const auto &intersection : intersections(basis.gridView(), e))
    {
      if (intersection.neighbor())
      {
        neighborLocalView.bind(intersection.outside());

        Dune::TypeTree::forEachLeafNode(
            localView.tree(),
            [&](const auto &insideNode, auto &&treePath)
            {
              const auto &outsideNode = Dune::TypeTree::child(neighborLocalView.tree(), treePath);

              std::vector<std::optional<int>> insideToOutside;
              insideToOutside.resize(insideNode.size());

              // Map all inside DOFs to outside DOFs if possible
              for (std::size_t i = 0; i < insideNode.size(); ++i)
              {
                for (std::size_t j = 0; j < outsideNode.size(); ++j)
                {
                  if (localView.index(insideNode.localIndex(i))
                      == neighborLocalView.index(outsideNode.localIndex(j)))
                  {
                    // Basis function should only appear once in the neighbor element.
                    test.check(not insideToOutside[i].has_value())
                        << "Basis function " << localView.index(insideNode.localIndex(i))
                        << " appears twice in element "
                        << elementStr(neighborLocalView.element(), basis.gridView());
                    insideToOutside[i] = j;
                  }
                }
              }

              // Apply continuity check on given intersection with given inside/outside DOF node
              // pair.
              auto isDifferentiable
                  = localCheck(intersection, treePath, insideNode, outsideNode, insideToOutside);

              for (std::size_t i = 0; i < insideNode.size(); ++i)
              {
                test.check(isDifferentiable[i])
                    << "Basis function " << localView.index(insideNode.localIndex(i))
                    << " is not differentiable " << flag.checkLocation
                    << " intersection of elements "
                    << elementStr(localView.element(), basis.gridView()) << " and "
                    << elementStr(neighborLocalView.element(), basis.gridView());
              }
            });
      }
    }
  }
  return test;
}
#endif
