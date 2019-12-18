// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/typetree/childextraction.hh>

#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/common/functionfromcallable.hh>
#include <dune/functions/common/functionconcepts.hh>

#include <dune/functions/backends/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/sizeinfo.hh>
#include <dune/functions/functionspacebases/flatvectorview.hh>
#include <dune/functions/functionspacebases/hierarchicnodetorangemap.hh>

#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>

namespace Dune {
namespace Functions {

namespace Imp {

struct AllTrueBitSetVector
{
  struct AllTrueBitSet
  {
    bool test(int i) const { return true; }
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



template <class B, class T, class NTRE, class HV, class LF, class HBV>
class LocalInterpolateVisitor
  : public TypeTree::TreeVisitor
  , public TypeTree::DynamicTraversal
{

public:

  using Basis = B;
  using LocalView = typename B::LocalView;
  using MultiIndex = typename LocalView::MultiIndex;

  using LocalFunction = LF;

  using Tree = T;

  using VectorBackend = HV;
  using BitVectorBackend = HBV;

  using NodeToRangeEntry = NTRE;

  using GridView = typename Basis::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using LocalDomain = typename Element::Geometry::LocalCoordinate;

  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  LocalInterpolateVisitor(const B& basis, HV& coeff, const HBV& bitVector, const LF& localF, const LocalView& localView, const NodeToRangeEntry& nodeToRangeEntry) :
    vector_(coeff),
    localF_(localF),
    bitVector_(bitVector),
    localView_(localView),
    nodeToRangeEntry_(nodeToRangeEntry)
  {
    static_assert(Dune::Functions::Concept::isCallable<LocalFunction, LocalDomain>(), "Function passed to LocalInterpolateVisitor does not model the Callable<LocalCoordinate> concept");
  }

  template<typename Node, typename TreePath>
  void pre(Node& node, TreePath treePath)
  {}

  template<typename Node, typename TreePath>
  void post(Node& node, TreePath treePath)
  {}

  template<typename Node, typename TreePath>
  void leaf(Node& node, TreePath treePath)
  {
    using FiniteElement = typename Node::FiniteElement;
    using FiniteElementRange = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
    using FiniteElementRangeField = typename FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType;
    using FunctionBaseClass = typename Dune::LocalFiniteElementFunctionBase<FiniteElement>::type;

    // Note that we capture j by reference. Hence we can switch
    // the selected component later on by modifying j. Maybe we
    // should avoid this naughty statefull lambda hack in favor
    // of a separate helper class.
    std::size_t j=0;
    auto localFj = [&](const LocalDomain& x){
      const auto& y = localF_(x);
      return flatVectorView(nodeToRangeEntry_(node, treePath, y))[j];
    };

    using FunctionFromCallable = typename Dune::Functions::FunctionFromCallable<FiniteElementRange(LocalDomain), decltype(localFj), FunctionBaseClass>;

    auto interpolationCoefficients = std::vector<FiniteElementRangeField>();

    auto&& fe = node.finiteElement();

    // We loop over j defined above and thus over the components of the
    // range type of localF_.

    auto blockSize = flatVectorView(vector_[localView_.index(0)]).size();

    for(j=0; j<blockSize; ++j)
    {

      // We cannot use localFj directly because interpolate requires a Dune::VirtualFunction like interface
      fe.localInterpolation().interpolate(FunctionFromCallable(localFj), interpolationCoefficients);
      for (size_t i=0; i<fe.localBasis().size(); ++i)
      {
        auto multiIndex = localView_.index(node.localIndex(i));
        auto bitVectorBlock = flatVectorView(bitVector_[multiIndex]);
        if (bitVectorBlock[j])
        {
          auto vectorBlock = flatVectorView(vector_[multiIndex]);
          vectorBlock[j] = interpolationCoefficients[i];
        }
      }
    }
  }


protected:

  VectorBackend& vector_;
  const LocalFunction& localF_;
  const BitVectorBackend& bitVector_;
  const LocalView& localView_;
  const NodeToRangeEntry& nodeToRangeEntry_;
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

  using Tree = typename B::LocalView::Tree;

  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  static_assert(Dune::Functions::Concept::isCallable<F, GlobalDomain>(), "Function passed to interpolate does not model the Callable<GlobalCoordinate> concept");

  auto&& gridView = basis.gridView();

  // Small helper functions to wrap vectors using istlVectorBackend
  // if they do not already satisfy the VectorBackend interface.
  auto toVectorBackend = [&](auto& v) -> decltype(auto) {
    return Dune::Hybrid::ifElse(models<Concept::VectorBackend<B>, decltype(v)>(),
    [&](auto id) -> decltype(auto) {
      return id(v);
    }, [&](auto id) -> decltype(auto) {
      return istlVectorBackend(id(v));
    });
  };
  auto toConstVectorBackend = [&](auto& v) -> decltype(auto) {
    return Dune::Hybrid::ifElse(models<Concept::ConstVectorBackend<B>, decltype(v)>(),
    [&](auto id) -> decltype(auto) {
      return id(v);
    }, [&](auto id) -> decltype(auto) {
      return istlVectorBackend(id(v));
    });
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

    Imp::LocalInterpolateVisitor<B, Tree, NTRE, decltype(vector), decltype(localF), decltype(bitVector)> localInterpolateVisitor(basis, vector, bitVector, localF, localView, nodeToRangeEntry);
    TypeTree::applyToTree(localView.tree(),localInterpolateVisitor);
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
