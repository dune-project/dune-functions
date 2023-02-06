// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/typetree/traversal.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/common/functionconcepts.hh>

#include <dune/functions/backends/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/sizeinfo.hh>
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



template<class VectorBackend, class BitVectorBackend, class LocalFunction, class LocalView, class NodeToRangeEntry>
void interpolateLocal(VectorBackend& vector, const BitVectorBackend& bitVector, const LocalFunction& localF, const LocalView& localView, const NodeToRangeEntry& nodeToRangeEntry)
{
  Dune::TypeTree::forEachLeafNode(localView.tree(), [&](auto&& node, auto&& treePath) {
    using Node = std::decay_t<decltype(node)>;
    using FiniteElement = typename Node::FiniteElement;
    using FiniteElementRange = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
    using FiniteElementRangeField = typename FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType;
    using LocalDomain = typename FiniteElement::Traits::LocalBasisType::Traits::DomainType;

    auto interpolationCoefficients = std::vector<FiniteElementRangeField>();
    auto&& fe = node.finiteElement();

    // for all other finite elements: use the FiniteElementRange directly for the interpolation
    auto localF_RE = [&](const LocalDomain& x){
      const auto& y = localF(x);
      return FiniteElementRange(nodeToRangeEntry(node, treePath, y));
    };

    fe.localInterpolation().interpolate(localF_RE, interpolationCoefficients);
    for (size_t i=0; i<fe.localBasis().size(); ++i)
    {
      auto multiIndex = localView.index(node.localIndex(i));
      if ( bitVector[multiIndex] )
        vector[multiIndex] = interpolationCoefficients[i];
    }
  });
}


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
  vector.resize(sizeInfo(basis));

  // Make a grid function supporting local evaluation out of f
  auto gf = makeGridViewFunction(f, gridView);

  // Obtain a local view of f
  auto localF = localFunction(gf);

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
