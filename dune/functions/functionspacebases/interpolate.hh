// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/typetree/childextraction.hh>

#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/common/functionfromcallable.hh>
#include <dune/functions/common/functionconcepts.hh>

#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/functionspacebases/sizeinfo.hh>
#include <dune/functions/functionspacebases/flatvectorbackend.hh>
#include <dune/functions/functionspacebases/defaultnodetorangemap.hh>

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
  using LocalIndexSet = typename B::LocalIndexSet;
  using MultiIndex = typename LocalIndexSet::MultiIndex;

  using LocalFunction = LF;

  using Tree = T;

  using HierarchicVector = HV;
  using HierarchicBitVector = HBV;

  using NodeToRangeEntry = NTRE;

  using GridView = typename Basis::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using LocalDomain = typename Element::Geometry::LocalCoordinate;

  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  using CoefficientBlock = typename std::decay<decltype(std::declval<HierarchicVector>()[std::declval<MultiIndex>()])>::type;
  using BitVectorBlock = typename std::decay<decltype(std::declval<HierarchicBitVector>()[std::declval<MultiIndex>()])>::type;

  LocalInterpolateVisitor(const B& basis, HV& coeff, const HBV& bitVector, const LF& localF, const LocalIndexSet& localIndexSet, const NodeToRangeEntry& nodeToRangeEntry) :
    vector_(coeff),
    localF_(localF),
    bitVector_(bitVector),
    localIndexSet_(localIndexSet),
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
    using FunctionBaseClass = typename Dune::LocalFiniteElementFunctionBase<FiniteElement>::type;

    // Note that we capture j by reference. Hence we can switch
    // the selected component later on by modifying j. Maybe we
    // should avoid this naughty statefull lambda hack in favor
    // of a separate helper class.
    std::size_t j=0;
    auto localFj = [&](const LocalDomain& x){
      const auto& y = localF_(x);
      const auto& y_node = nodeToRangeEntry_(node, y);
      using FunctionRange = typename std::decay<decltype(y_node)>::type;
      return FlatVectorBackend<FunctionRange>::getEntry(y_node, j);
    };

    using FunctionFromCallable = typename Dune::Functions::FunctionFromCallable<FiniteElementRange(LocalDomain), decltype(localFj), FunctionBaseClass>;

    auto interpolationValues = std::vector<FiniteElementRange>();

    auto&& fe = node.finiteElement();

    // We loop over j defined above and thus over the components of the
    // range type of localF_.

    auto blockSize = FlatVectorBackend<CoefficientBlock>::size(vector_[localIndexSet_.index(0)]);

    for(j=0; j<blockSize; ++j)
    {

      // We cannot use localFj directly because interpolate requires a Dune::VirtualFunction like interface
      fe.localInterpolation().interpolate(FunctionFromCallable(localFj), interpolationValues);
      for (size_t i=0; i<fe.localBasis().size(); ++i)
      {
        auto multiIndex = localIndexSet_.index(node.localIndex(i));
        const auto& bitVectorBlock = bitVector_[multiIndex];
        const auto& interpolateHere = FlatVectorBackend<BitVectorBlock>::getEntry(bitVectorBlock,j);

        if (interpolateHere)
        {
          auto&& vectorBlock = vector_[multiIndex];
          FlatVectorBackend<CoefficientBlock>::getEntry(vectorBlock, j) = interpolationValues[i];
        }
      }
    }
  }


protected:

  HierarchicVector& vector_;
  const LocalFunction& localF_;
  const HierarchicBitVector& bitVector_;
  const LocalIndexSet& localIndexSet_;
  const NodeToRangeEntry& nodeToRangeEntry_;
};


} // namespace Imp




/**
 * \brief Interpolate given function in discrete function space
 *
 * Interpolation is done wrt the leaf node of the ansatz tree
 * corresponding to the given tree path.
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * FlatVectorBackend.
 *
 * \param basis Global function space basis of discrete function space
 * \param treePath Tree path specifying the part of the ansatz tree to use
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 * \param nodeToRangeEntry Polymorphic functor mapping local ansatz nodes to range-indices of given function
 * \param bitVector A vector with flags marking all DOFs that should be interpolated
 */
template <class B, class... TreeIndices, class NTRE, class C, class F, class BV>
void interpolateTreeSubset(const B& basis, const TypeTree::HybridTreePath<TreeIndices...>& treePath, C&& coeff, const F& f, const NTRE& nodeToRangeEntry, const BV& bv)
{
  using GridView = typename B::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using Tree = typename std::decay<decltype(TypeTree::child(basis.localView().tree(),treePath))>::type;

  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  static_assert(Dune::Functions::Concept::isCallable<F, GlobalDomain>(), "Function passed to interpolate does not model the Callable<GlobalCoordinate> concept");

  auto&& gridView = basis.gridView();

  auto&& bitVector = makeHierarchicVectorForMultiIndex<typename B::MultiIndex>(bv);
  auto&& vector = makeHierarchicVectorForMultiIndex<typename B::MultiIndex>(coeff);
  vector.resize(sizeInfo(basis));



  // Make a grid function supporting local evaluation out of f
  auto gf = makeGridViewFunction(f, gridView);

  // Obtain a local view of f
  auto localF = localFunction(gf);

  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();

  for (const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);
    localF.bind(e);

    auto&& subTree = TypeTree::child(localView.tree(),treePath);

    Imp::LocalInterpolateVisitor<B, Tree, NTRE, decltype(vector), decltype(localF), decltype(bitVector)> localInterpolateVisitor(basis, vector, bitVector, localF, localIndexSet, nodeToRangeEntry);
    TypeTree::applyToTree(subTree,localInterpolateVisitor);
  }
}


template <class B, class... TreeIndices, class C, class F, class BV>
void interpolateTreeSubset(const B& basis, const TypeTree::HybridTreePath<TreeIndices...>& treePath, C&& coeff, const F& f, const BV& bitVector)
{
  interpolateTreeSubset(basis, treePath, coeff, f, makeDefaultNodeToRangeMap(basis, treePath), bitVector);
}


template <class B, class... TreeIndices, class NTRE, class C, class F>
void interpolateTree(const B& basis, const TypeTree::HybridTreePath<TreeIndices...>& treePath, C&& coeff, const F& f, const NTRE& nodeToRangeEntry)
{
  interpolateTreeSubset(basis, treePath, coeff, f, nodeToRangeEntry, Imp::AllTrueBitSetVector());
}


template <class B, class... TreeIndices, class C, class F>
void interpolateTree(const B& basis, const TypeTree::HybridTreePath<TreeIndices...>& treePath, C&& coeff, const F& f)
{
  interpolateTreeSubset(basis, treePath, coeff, f, makeDefaultNodeToRangeMap(basis, treePath), Imp::AllTrueBitSetVector());
}



/**
 * \brief Interpolate given function in discrete function space
 *
 * Interpolation is done wrt the leaf node of the ansatz tree
 * corresponding to the given tree path.
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * FlatVectorBackend.
 *
 * \param basis Global function space basis of discrete function space
 * \param treePath Tree path specifying the part of the ansatz tree to use
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 * \param bitVector A vector with flags marking all DOFs that should be interpolated
 */
template <class B, class... TreeIndices, class C, class F, class BV>
void interpolate(const B& basis, const TypeTree::HybridTreePath<TreeIndices...>& treePath, C&& coeff, const F& f, const BV& bitVector)
{
  interpolateTreeSubset(basis, treePath, coeff, f, makeDefaultNodeToRangeMap(basis, treePath), bitVector);
}

namespace Imp {

  template<class... T>
  constexpr std::true_type isHybridTreePath(const TypeTree::HybridTreePath<T...>&)
  { return {}; }

  template<class T>
  constexpr std::false_type isHybridTreePath(const T&)
  { return {}; }

  template<class T>
  constexpr auto isHybridTreePath() -> decltype(isHybridTreePath(std::declval<std::decay_t<T>>()))
  { return {}; }

}

/**
 * \brief Interpolate given function in discrete function space
 *
 * Interpolation is done wrt the leaf node of the ansatz tree
 * corresponding to the given tree path.  Only vector coefficients marked as 'true' in the
 * bitVector argument are interpolated.  Use this, e.g., to interpolate Dirichlet boundary values.
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * FlatVectorBackend.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 * \param bitVector A vector with flags marking all DOFs that should be interpolated
 */
template <class B, class C, class F, class BV,
         std::enable_if_t<not Imp::isHybridTreePath<C>(), int> = 0>
void interpolate(const B& basis, C&& coeff, const F& f, const BV& bitVector)
{
  auto root = Dune::TypeTree::hybridTreePath();
  interpolateTreeSubset(basis, root, coeff, f, makeDefaultNodeToRangeMap(basis, root), bitVector);
}



/**
 * \brief Interpolate given function in discrete function space
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * FlatVectorBackend.
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
  interpolate (basis, Dune::TypeTree::hybridTreePath(), coeff, f, Imp::AllTrueBitSetVector());
}

/**
 * \brief Interpolate given function in discrete function space
 *
 * Interpolation is done wrt the leaf node of the ansatz tree
 * corresponding to the given tree path.
 *
 * Notice that this will only work if the range type of f and
 * the block type of corresponding coeff entries are compatible
 * and supported by FlatVectorBackend.
 *
 * \param basis Global function space basis of discrete function space
 * \param treePath Tree path specifying the part of the ansatz tree to use
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 */
template <class B, class... TreeIndices, class C, class F>
void interpolate(const B& basis, const TypeTree::HybridTreePath<TreeIndices...>& treePath, C&& coeff, const F& f)
{
  interpolate (basis, treePath, coeff, f, Imp::AllTrueBitSetVector());
}

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
