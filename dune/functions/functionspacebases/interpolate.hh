// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/referencehelper.hh>

#include <dune/typetree/traversal.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/common/functionconcepts.hh>

#include <dune/functions/backends/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/flatvectorview.hh>
#include <dune/functions/functionspacebases/hierarchicnodetorangemap.hh>

namespace Dune {
namespace Functions {

namespace Imp {

struct AllTrueBitSetVector
{
  struct AllTrueBitSet
  {
    bool test(int) const { return true; }
  } allTrue_;

  operator bool() const
  {
    return true;
  }

  template<class I>
  const AllTrueBitSetVector& operator[](const I&) const
  {
    return *this;
  }

  template<class SP>
  void resize(const SP&) const
  {}

};



// This helper function implements the restriction of some given function of type F.
// The restriction is a simple callback that is applied to the values of the
// function and the values of its derivative.
template<class F, class Restriction>
class ComponentFunction
{
public:

  ComponentFunction(F f, Restriction restriction) :
    f_(std::move(f)),
    restriction_(std::move(restriction))
  {}

  template<class Domain>
  auto operator()(const Domain& x) const
  {
    return restriction_(f_(x));
  }

  friend auto derivative(const ComponentFunction& cf)
  {
    // This provides support for capturing the derivative of the function by reference
    // using forwardCapture for perfect forwarding capture. If the function caches its
    // derivative, this saves a potentially costly copy.
    auto&& df = derivative(Dune::resolveRef(cf.f_));
    return [&, df=forwardCapture(std::forward<decltype(df)>(df))](auto&& x) {
      return cf.restriction_(df.forward()(x));
    };
  }

private:
  F f_;
  Restriction restriction_;
};




// This helper function implements caching of the derivative for local functions.
// When using an algorithm that gets a LocalFunction and calls its derivative
// on each element, this leads to a costly call of derivative(f). E.g. for a
// DiscreteGlobalBasisFunction, this will allocate several buffer.
// To avoid this, this helper function caches the derivative and hands
// out the cached derivative by reference. To ensure that the handed out
// derivative is appropriately bound, binding the function will automatically
// bind the cached derivative.
//
// Notice that we cannot simply create the derivative in the constructor,
// because this may throw for functions that do not implement the derivative.
template<class F>
class CachedDerivativeLocalFunction
{
  using Derivative = std::decay_t<decltype(derivative(Dune::resolveRef(std::declval<const F&>())))>;

public:

  CachedDerivativeLocalFunction(F f) :
    f_(f)
  {}

  template<class Element>
  void bind(const Element& element)
  {
    Dune::resolveRef(f_).bind(element);
    if (derivative_)
      derivative_.value().bind(element);
  }

  template<class X>
  auto operator()(const X& x) const
  {
    return f_(x);
  }

  friend const Derivative& derivative(const CachedDerivativeLocalFunction& cdlf)
  {
    if (not cdlf.derivative_)
    {
      auto&& lf = Dune::resolveRef(cdlf.f_);
      cdlf.derivative_ = derivative(lf);
      if (lf.bound())
        cdlf.derivative_.value().bind(lf.localContext());
    }
    return cdlf.derivative_.value();
  }

private:
  F f_;
  mutable std::optional<Derivative> derivative_;
};



template<class VectorBackend, class BitVectorBackend, class LocalFunction, class LocalView, class NodeToRangeEntry>
void interpolateLocal(VectorBackend& vector, const BitVectorBackend& bitVector, const LocalFunction& localF, const LocalView& localView, const NodeToRangeEntry& nodeToRangeEntry)
{
  Dune::TypeTree::forEachLeafNode(localView.tree(), [&](auto&& node, auto&& treePath) {
    if (node.empty())
      return;
    using Node = std::decay_t<decltype(node)>;
    using FiniteElement = typename Node::FiniteElement;
    using FiniteElementRangeField = typename FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType;

    auto interpolationCoefficients = std::vector<FiniteElementRangeField>();
    auto&& fe = node.finiteElement();
    auto localF_RE = ComponentFunction(std::cref(localF), [&](auto&& y) { return nodeToRangeEntry(node, treePath, y); });

    fe.localInterpolation().interpolate(localF_RE, interpolationCoefficients);
    for (size_t i=0; i<fe.localBasis().size(); ++i)
    {
      auto multiIndex = localView.index(node.localIndex(i));
      if ( bitVector[multiIndex] )
        vector[multiIndex] = interpolationCoefficients[i];
    }
  });
}


struct HasDerivative
{
  template<class F>
  auto require(F&& f) -> decltype(derivative(f));
};

} // namespace Imp




/**
 * \brief Interpolate given function in discrete function space
 *
 * Only vector coefficients marked as 'true' in the
 * bitVector argument are interpolated.  Use this, e.g., to
 * interpolate Dirichlet boundary values.
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * flatVectorView.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 * \param bitVector A vector with flags marking all DOFs that should be interpolated
 * \param nodeToRangeEntry Polymorphic functor mapping local ansatz nodes to range-indices of given function
 */
template <class B, class C, class F, class BV, class NTRE>
void interpolate(const B& basis, C&& coeff, const F& f, const BV& bv, const NTRE& nodeToRangeEntry)
{
  using GridView = typename B::GridView;
  using Element = typename GridView::template Codim<0>::Entity;
  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  static_assert(Dune::Functions::Concept::isCallable<F, GlobalDomain>(), "Function passed to interpolate does not model the Callable<GlobalCoordinate> concept");

  auto&& gridView = basis.gridView();

  // Small helper functions to wrap vectors using istlVectorBackend
  // if they do not already satisfy the VectorBackend interface.
  auto toVectorBackend = [&](auto& v) -> decltype(auto) {
    if constexpr (models<Concept::VectorBackend<B>, decltype(v)>()) {
      return v;
    } else {
      return istlVectorBackend(v);
    }
  };

  auto toConstVectorBackend = [&](auto& v) -> decltype(auto) {
    if constexpr (models<Concept::ConstVectorBackend<B>, decltype(v)>()) {
      return v;
    } else {
      return istlVectorBackend(v);
    }
  };

  auto&& bitVector = toConstVectorBackend(bv);
  auto&& vector = toVectorBackend(coeff);
  vector.resize(basis);

  // Make a grid function supporting local evaluation out of f
  auto gf = makeGridViewFunction(f, gridView);

  // Obtain a local view of f
  // To avoid costly reconstruction of the derivative on each element,
  // we use the CachedDerivativeLocalFunction wrapper if the function
  // is differentiable. This wrapper will handout
  // a reference to a single cached derivative object.
  auto localF = [&](){
    if constexpr (models<Imp::HasDerivative, decltype(localFunction(gf))>())
      return Imp::CachedDerivativeLocalFunction(localFunction(gf));
    else
      return localFunction(gf);
  }();

  auto localView = basis.localView();

  for (const auto& e : elements(gridView))
  {
    localView.bind(e);
    localF.bind(e);
    Imp::interpolateLocal(vector, bitVector, localF, localView, nodeToRangeEntry);
  }
}

/**
 * \brief Interpolate given function in discrete function space
 *
 * Only vector coefficients marked as 'true' in the
 * bitVector argument are interpolated.  Use this, e.g., to
 * interpolate Dirichlet boundary values.
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * flatVectorView.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 * \param bitVector A vector with flags marking all DOFs that should be interpolated
 */
template <class B, class C, class F, class BV>
void interpolate(const B& basis, C&& coeff, const F& f, const BV& bitVector)
{
  interpolate(basis, coeff, f, bitVector, HierarchicNodeToRangeMap());
}

/**
 * \brief Interpolate given function in discrete function space
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * flatVectorView.
 *
 * This function will only work, if the local ansatz tree of
 * the basis is trivial, i.e., a single leaf node.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 */
template <class B, class C, class F>
void interpolate(const B& basis, C&& coeff, const F& f)
{
  interpolate (basis, coeff, f, Imp::AllTrueBitSetVector(), HierarchicNodeToRangeMap());
}

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
