// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/common/functionfromcallable.hh>
#include <dune/functions/common/functionconcepts.hh>

#include <dune/functions/functionspacebases/hierarchicvectorbackend.hh>
#include <dune/functions/functionspacebases/flatvectorbackend.hh>

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

  const AllTrueBitSetVector& operator[](int i) const
  {
    return *this;
  }

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
 * \param bitVector A vector with flags marking ald DOFs that should be interpolated
 */
template <class B, class TP, class C, class F, class BV>
void interpolate(const B& basis, TP&& treePath, C& coeff, F&& f, BV&& bitVector)
{
  using GridView = typename B::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using Tree = typename std::decay<decltype(getChild(basis.localView().tree(), treePath))>::type;
  using FiniteElement = typename Tree::FiniteElement;
  using FunctionBaseClass = typename Dune::LocalFiniteElementFunctionBase<FiniteElement>::type;

  using LocalBasisRange = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
  using LocalDomain = typename Element::Geometry::LocalCoordinate;

  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  using Backend = Dune::Functions::HierarchicVectorBackend;

  using CoefficientBlock = typename std::decay<decltype(
                                                        Backend::getEntry(coeff, basis.indexSet().localIndexSet().index(0))
                                                        )>::type;
  using BitVectorBlock = typename std::decay<decltype(
                                                        Backend::getEntry(bitVector, basis.indexSet().localIndexSet().index(0))
                                                        )>::type;

  static_assert(Dune::Functions::Concept::isCallable<F, GlobalDomain>(), "Function passed to interpolate does not model the Callable<GlobalCoordinate> concept");

  auto&& gridView = basis.gridView();

  // Make a grid function supporting local evaluation out of f
  auto gf = makeGridViewFunction(std::forward<F>(f), gridView);

  // Obtain a local view of f
  auto localF = localFunction(gf);

  // Note that we capture j by reference. Hence we can switch
  // the selected component later on by modifying j. Maybe we
  // should avoid this naughty statefull lambda hack in favor
  // of a separate helper class.
  int j=0;
  auto localFj = [&](const LocalDomain& x){
    using FunctionRange = typename std::decay<decltype(localF(LocalDomain(0)))>::type;
    return FlatVectorBackend<FunctionRange>::getEntry(localF(x), j);
  };

  using FunctionFromCallable = typename Dune::Functions::FunctionFromCallable<LocalBasisRange(LocalDomain), decltype(localFj), FunctionBaseClass>;

  auto basisIndexSet = basis.indexSet();
  coeff.resize(basisIndexSet.size());

//  auto processed = Dune::BitSetVector<1>(basisIndexSet.size(), false);
  auto interpolationValues = std::vector<LocalBasisRange>();

  auto localView = basis.localView();
  auto localIndexSet = basisIndexSet.localIndexSet();

//  auto blockSize = Imp::FlatIndexContainerAccess<CoefficientBlock>::size(coeff[0]);

  for (const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);
    localF.bind(e);

    auto&& node = getChild(localView.tree(), treePath);
    auto&& fe = node.finiteElement();

#if 0
    // check if all components have already been processed
    bool allProcessed = true;
    for (size_t i=0; i<fe.localBasis().size(); ++i)
    {
      // if index was already processed we don't need to do further checks
      auto index = localIndexSet.index(i)[0];
      if (processed[index][0])
        continue;

      // if index was not processed, check if any entry is marked for interpolation
      auto&& bitVectorBlock = bitVector[index];
      for(int k=0; k<blockSize; ++k)
      {
        if (Imp::FlatIndexContainerAccess<BitVectorBlock>::getEntry(bitVectorBlock,k))
        {
          allProcessed = false;
          break;
        }
      }
    }

    if (not(allProcessed))
#endif
    {
      // We loop over j defined above and thus over the components of the
      // range type of localF.

      auto blockSize = FlatVectorBackend<CoefficientBlock>::size(Backend::getEntry(coeff, localIndexSet.index(0)));

      for(j=0; j<blockSize; ++j)
      {

        // We cannot use localFj directly because interpolate requires a Dune::VirtualFunction like interface
        fe.localInterpolation().interpolate(FunctionFromCallable(localFj), interpolationValues);
        for (size_t i=0; i<fe.localBasis().size(); ++i)
        {
          auto multiIndex = localIndexSet.index(node.localIndex(i));
          const auto& bitVectorBlock = Backend::getEntry(bitVector, multiIndex);
          const auto& interpolateHere = FlatVectorBackend<BitVectorBlock>::getEntry(bitVectorBlock,j);

//          if (not(processed[index][0]) and interpolateHere)
          if (interpolateHere)
          {
            auto&& vectorBlock = Backend::getEntry(coeff, multiIndex);
            FlatVectorBackend<CoefficientBlock>::getEntry(vectorBlock, j) = interpolationValues[i];
          }
        }
      }
//      for (size_t i=0; i<fe.localBasis().size(); ++i)
//        processed[localIndexSet.index(i)[0]][0] = true;
    }
  }
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
void interpolate(const B& basis, C& coeff, F&& f)
{
  interpolate (basis, std::make_tuple(), coeff, f, Imp::AllTrueBitSetVector());
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
template <class B, class TreePath, class C, class F>
void interpolate(const B& basis, TreePath&& treePath, C& coeff, F&& f)
{
  interpolate (basis, treePath, coeff, f, Imp::AllTrueBitSetVector());
}

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
