// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CUBICHERMITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CUBICHERMITEBASIS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/tuplevector.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localfiniteelementvariant.hh>
#include <dune/localfunctions/rannacherturek.hh>
#include <dune/localfunctions/crouzeixraviart.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

namespace Impl {

  // Helper function returning a random access range
  // of global indices associated to the element. The ordering
  // is derived from the LocalCoefficients object.
  // Having this as a member of MCMGMapper would be nice.
  // But this would introduce the LocalCoefficient in dune-grid.
  // This is at least a soft, 'conceptual' dependency issue.
  template<class GridView, class LocalCoefficients>
  auto subIndexRange(const Dune::MultipleCodimMultipleGeomTypeMapper<GridView>& mapper, const typename GridView::template Codim<0>::Entity& element, const LocalCoefficients& localCoefficients)
  {
    // Here we make use of the 'hidden' (poorly documented) MCMGMapper feature to support
    // multiple DOFs per subentity. However, we do not take care for any reordering.
    return Dune::transformedRangeView(Dune::range(localCoefficients.size()), [&](auto localIndex) {
      auto localKey = localCoefficients.localKey(localIndex);
      return mapper.subIndex(element, localKey.subEntity(), localKey.codim()) + localKey.index();
    });
  }

  // Helper function returning an unordered range
  // of global indices associated to the element.
  // This could be implemented cheaper internally in
  // the MCMGMapper by storing a precomputed
  // container of all subsentities addressed by the layout.
  template<class GridView>
  auto subIndexSet(const Dune::MultipleCodimMultipleGeomTypeMapper<GridView>& mapper, const typename GridView::template Codim<0>::Entity& element)
  {
    using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using Index = typename Mapper::Index;
    constexpr auto dimension = GridView::dimension;
    auto subIndices = std::vector<Index>();
    auto referenceElement = Dune::referenceElement<double, dimension>(element.type());
    for(auto codim : Dune::range(dimension+1))
    {
      for(auto subEntity : Dune::range(referenceElement.size(codim)))
      {
        std::size_t c = mapper.layout()(referenceElement.type(subEntity, codim), dimension);
        if (c>0)
        {
          std::size_t firstIndex = mapper.subIndex(element, subEntity, codim);
          for(auto j : Dune::range(firstIndex, firstIndex+c))
          {
            subIndices.push_back(j);
          }
        }
      }
    }
    return subIndices;
  }

  // Helper function computing an average mesh size per subentity
  // by averaging over the adjacent elements. This only considers
  // the subentities handled by the given mapper and returns a
  // vector of mesh sizes indixed according to the mapper.
  template<class Mapper>
  auto computeAverageSubEntityMeshSize(const Mapper& mapper)
  {
    constexpr auto dimension = Mapper::GridView::dimension;

    std::vector<unsigned int> adjacentElements(mapper.size(), 0);
    std::vector<double> subEntityMeshSize(mapper.size(), 0.0);
    for(const auto& element : Dune::elements(mapper.gridView()))
    {
      auto A = element.geometry().volume();
      for(auto i : Impl::subIndexSet(mapper, element))
      {
        subEntityMeshSize[i] += A;
        ++(adjacentElements[i]);
      }
    }
    for(auto i : Dune::range(mapper.size()))
      subEntityMeshSize[i] = std::pow(subEntityMeshSize[i]/adjacentElements[i], 1./dimension);
    return subEntityMeshSize;
  }



} // namespace Impl



/**
 * \brief A generic MixIn class for PreBasis with flat indices computed from a mapper.
 *
 * This abstracts all index computations that can be implemented using a
 * MultipleCodimMultipleGeomTypeMapper with appropriate MCMGLayout..
 * In order to use this, you need to derive from this class and
 * pass the layout in the constructor. Then the mixin takes care
 * for all the index and size computation and the derived class
 * only needs to add the node creation.
 *
 * Be careful: This does not do any reordering of indices
 * if multiple basis functions are associated to the same
 * subentity.
 *
 * \tparam GV The grid view the basis is defined on.
 */
template<typename GV>
class LeafPreBasisMapperMixIn
{
  static const int gridDim = GV::dimension;

public:

  using GridView = GV;
  using size_type = std::size_t;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  LeafPreBasisMapperMixIn(const GridView& gv, Dune::MCMGLayout layout) :
    gridView_(gv),
    mapper_(gridView_, std::move(layout))
  {}

  void initializeIndices()
  {
    maxNodeSize_ = 0;
    for(const GeometryType& elementType : gridView_.indexSet().types(0))
    {
      auto referenceElement = Dune::referenceElement<double, gridDim>(elementType);
      for(auto codim : Dune::range(gridDim+1))
        for(auto i : Dune::range(referenceElement.size(codim)))
          maxNodeSize_ += mapper_.layout()(referenceElement.type(i, codim), gridDim);
    }
  }

  const GridView& gridView() const
  {
    return gridView_;
  }

  void update(const GridView& gv)
  {
    gridView_ = gv;
    mapper_.update(gridView_);
  }

  size_type size() const
  {
    return mapper_.size();
  }

  template<class SizePrefix>
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return maxNodeSize_;
  }

  template<class Node, class It>
  It indices(const Node& node, It it) const
  {
    for(const auto& globalIndex : Impl::subIndexRange(mapper_, node.element(), node.finiteElement().localCoefficients()))
    {
      *it = {{ (size_type)globalIndex }};
      ++it;
    }
    return it;
  }

protected:
  GridView gridView_;
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper_;
  std::size_t maxNodeSize_;
};



namespace Impl {


// *****************************************************************************
// * Some helper functions for building polynomial bases from monomials
// *****************************************************************************

// Evaluation of 1d monomial values
template<class K>
static constexpr auto evaluateMonomialValues(const Dune::FieldVector<K,1>& x)
{
  using Range = Dune::FieldVector<K,1>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t size = (maxOrder+1);
  auto xPowers = std::array<double,maxOrder+1>{};
  xPowers[0] = 1.0;
  for(auto k: Dune::range(maxOrder))
    xPowers[k+1] = xPowers[k]*x[0];
  auto y = Dune::FieldVector<Range,size>{};
  for(auto order : Dune::range(maxOrder+1))
    y[order] = xPowers[order];
  return y;
}

// Evaluation of 1d monomial jacobians
template<class K>
static constexpr auto evaluateMonomialJacobians(const Dune::FieldVector<K,1>& x)
{
  using Jacobian = Dune::FieldMatrix<K,1,1>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t size = (maxOrder+1);
  auto xPowers = std::array<double,maxOrder+1>{};
  xPowers[0] = 1.0;
  for(auto k: Dune::range(maxOrder))
    xPowers[k+1] = xPowers[k]*x[0];
  auto y = Dune::FieldVector<Jacobian,size>{};
  for(auto order : Dune::range(std::size_t(1), maxOrder+1))
    y[order][0][2] = order*xPowers[order-1];
  return y;
}

// Evaluation of 2d monomial values
template<class K>
static constexpr auto evaluateMonomialValues(const Dune::FieldVector<K,2>& x)
{
  using Range = Dune::FieldVector<K,1>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t dim=2;
  constexpr std::size_t size = (maxOrder+1)*(maxOrder+2)/2;
  auto xPowers = std::array<std::array<double,maxOrder+1>,dim>{};
  for(auto j: Dune::range(dim))
  {
    xPowers[j][0] = 1.0;
    for(auto k: Dune::range(maxOrder))
      xPowers[j][k+1] = xPowers[j][k]*x[j];
  }
  auto y = Dune::FieldVector<Range,size>{};
  std::size_t index=0;
  for(auto order : Dune::range(maxOrder+1))
  {
    for(auto k : Dune::range(order+1))
    {
      y[index] = xPowers[0][order-k]*xPowers[1][k];
      ++index;
    }
  }
  return y;
}

// Evaluation of 2d monomial jacobians
template<class K>
static constexpr auto evaluateMonomialJacobians(const Dune::FieldVector<K,2>& x)
{
  using Jacobian = Dune::FieldMatrix<K,1,2>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t dim=2;
  constexpr std::size_t size = (maxOrder+1)*(maxOrder+2)/2;
  auto xPowers = std::array<std::array<double,maxOrder+1>,dim>{};
  for(auto j: Dune::range(dim))
  {
    xPowers[j][0] = 1.0;
    for(auto k: Dune::range(maxOrder))
      xPowers[j][k+1] = xPowers[j][k]*x[j];
  }
  auto y = Dune::FieldVector<Jacobian,size>{};
  std::size_t index=0;
  for(auto order : Dune::range(maxOrder+1))
  {
    for(auto k : Dune::range(order+1))
    {
      if (order-k>0)
        y[index][0][0] = (order-k)*xPowers[0][order-k-1]*xPowers[1][k];
      if (k>0)
        y[index][0][1] = k*xPowers[0][order-k]*xPowers[1][k-1];
      ++index;
    }
  }
  return y;
}



// *****************************************************************************
// * CubicHermiteLocalFiniteElement
// *****************************************************************************



template<class DF, class RF, unsigned int dim, bool reduced>
class CubicHermiteLocalBasis
{

  static constexpr auto makeReferenceBasisCoefficients() {
    if constexpr (dim==1)
      return Dune::FieldMatrix<int,4,4>{
        { 1,    0,   -3,    2},
        { 0,    1,   -2,    1},
        { 0,    0,    3,   -2},
        { 0,    0,   -1,    1}
      };
    if constexpr ((dim==2) and (not reduced))
      return Dune::FieldMatrix<int,10,10>{
        { 1,    0,    0,   -3,  -13,   -3,    2,   13,   13,    2},
        { 0,    1,    0,   -2,   -3,    0,    1,    3,    2,    0},
        { 0,    0,    1,    0,   -3,   -2,    0,    2,    3,    1},
        { 0,    0,    0,    3,   -7,    0,   -2,    7,    7,    0},
        { 0,    0,    0,   -1,    2,    0,    1,   -2,   -2,    0},
        { 0,    0,    0,    0,   -1,    0,    0,    2,    1,    0},
        { 0,    0,    0,    0,   -7,    3,    0,    7,    7,   -2},
        { 0,    0,    0,    0,   -1,    0,    0,    1,    2,    0},
        { 0,    0,    0,    0,    2,   -1,    0,   -2,   -2,    1},
        { 0,    0,    0,    0,   27,    0,    0,  -27,  -27,    0}
      };
    if constexpr ((dim==2) and (reduced))
    {
      auto w = std::array{1./3, 1./18, 1./18, 1./3, -1./9, 1./18, 1./3, 1./18, -1./9};
      return Dune::FieldMatrix<double,9,10>{
        { 1,    0,    0,   -3,  -13 + w[0]*27,   -3,    2,   13 - w[0]*27,   13 - w[0]*27,    2},
        { 0,    1,    0,   -2,   -3 + w[1]*27,    0,    1,    3 - w[1]*27,    2 - w[1]*27,    0},
        { 0,    0,    1,    0,   -3 + w[2]*27,   -2,    0,    2 - w[2]*27,    3 - w[2]*27,    1},
        { 0,    0,    0,    3,   -7 + w[3]*27,    0,   -2,    7 - w[3]*27,    7 - w[3]*27,    0},
        { 0,    0,    0,   -1,    2 + w[4]*27,    0,    1,   -2 - w[4]*27,   -2 - w[4]*27,    0},
        { 0,    0,    0,    0,   -1 + w[5]*27,    0,    0,    2 - w[5]*27,    1 - w[5]*27,    0},
        { 0,    0,    0,    0,   -7 + w[6]*27,    3,    0,    7 - w[6]*27,    7 - w[6]*27,   -2},
        { 0,    0,    0,    0,   -1 + w[7]*27,    0,    0,    1 - w[7]*27,    2 - w[7]*27,    0},
        { 0,    0,    0,    0,    2 + w[8]*27,   -1,    0,   -2 - w[8]*27,   -2 - w[8]*27,    1},
      };
    }
  }

  // These are the coefficients of the cubic Hermite basis functions
  // on the reference element wrt. the monomials. These have been computed
  // by solving with the corresponiding Vandermonde-matrix for the reference
  // element in advance.
  // The basis functions can be evaluated by first evaluating the monomials
  // and then transforming their values with these coefficients using
  //
  //   referenceBasisCoefficients.mv(evaluateMonomialValues(x), values);
  //
  // Surprisingly, storing them as static constexpr member is slightly faster
  // than using the static constexpr function directly.
  static constexpr auto referenceBasisCoefficients = makeReferenceBasisCoefficients();

  // This transforms the function or derivative values from the basis
  // functions on the reference element to those on the grid element.
  // Since Hermite elements do not form an affine family, the transformation
  // of derivative DOFs involves the Jacobian of the grid element transformation.
  // To avoid blowup of condition numbers due to h-dependend basis functions,
  // an h-dependent rescaling is applied.
  template<class LambdaRefValues, class Entry>
  void transformToElementBasis(const LambdaRefValues& refValues, std::vector<Entry>& out) const
  {
    if constexpr (dim==1)
    {
      const auto& J = elementJacobian_;
      out.resize(refValues.size());
      out[0] = refValues[0];
      out[1] = J*refValues[1] / (*localSubEntityMeshSize_)[1];
      out[2] = refValues[2];
      out[3] = J*refValues[3] / (*localSubEntityMeshSize_)[1];;
    }
    if constexpr (dim==2)
    {
      const auto& J = elementJacobian_;
      out.resize(refValues.size());
      out[0] = refValues[0];
      out[1] = (J[0][0]*refValues[1] + J[0][1]*refValues[2]) / (*localSubEntityMeshSize_)[1];
      out[2] = (J[1][0]*refValues[1] + J[1][1]*refValues[2]) / (*localSubEntityMeshSize_)[2];
      out[3] = refValues[3];
      out[4] = (J[0][0]*refValues[4] + J[0][1]*refValues[5]) / (*localSubEntityMeshSize_)[4];
      out[5] = (J[1][0]*refValues[4] + J[1][1]*refValues[5]) / (*localSubEntityMeshSize_)[5];
      out[6] = refValues[6];
      out[7] = (J[0][0]*refValues[7] + J[0][1]*refValues[8]) / (*localSubEntityMeshSize_)[7];
      out[8] = (J[1][0]*refValues[7] + J[1][1]*refValues[8]) / (*localSubEntityMeshSize_)[8];
      if constexpr (not reduced)
        out[9] = refValues[9];
    }
  }

  using ElementJacobian = Dune::FieldMatrix<DF, dim,dim>;

  using LocalSubEntityMeshSize = std::vector<double>;

public:

  using Domain = Dune::FieldVector<DF, dim>;
  using Range = Dune::FieldVector<RF, 1>;
  using Jacobian = Dune::FieldMatrix<RF, 1, dim>;
  using Traits = Dune::LocalBasisTraits<DF, dim, Domain, RF, 1, Range, Jacobian>;
  using OrderArray = std::array<unsigned int, dim>;

  static constexpr unsigned int size()
  {
    return decltype(referenceBasisCoefficients)::rows;
  }

  inline void evaluateFunction(const Domain& x, std::vector<Range>& values) const
  {
    auto monomialValues = evaluateMonomialValues(x);
    auto referenceValues = Dune::FieldVector<Range, size()>{};
    referenceBasisCoefficients.mv(monomialValues, referenceValues);
    transformToElementBasis(referenceValues, values);
  }

  inline void evaluateJacobian(const Domain& x, std::vector<Jacobian>& jacobians) const
  {
    auto monomialJacobians = evaluateMonomialJacobians(x);
    auto referenceJacobians = Dune::FieldVector<Jacobian, size()>{};
    referenceBasisCoefficients.mv(monomialJacobians, referenceJacobians);
    transformToElementBasis(referenceJacobians, jacobians);
  }

  void partial(const OrderArray& order, const Domain& x, std::vector<Range>& out) const
  {
    auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
    if (totalOrder == 0)
      evaluateFunction(x, out);
    DUNE_THROW(RangeError, "partial() not implemented for given order");
  }

  unsigned int order() const
  {
    return 3;
  }

  template<class Element>
  void bind(const Element& element, const LocalSubEntityMeshSize& localSubEntityMeshSize) {
    localSubEntityMeshSize_ = &localSubEntityMeshSize;
    auto center = Dune::ReferenceElements<DF, dim>::simplex().position(0, 0);
    elementJacobian_ = element.geometry().jacobian(center);
  }

private:
  ElementJacobian elementJacobian_;
  const LocalSubEntityMeshSize* localSubEntityMeshSize_;
};



template<unsigned int dim, bool reduced>
struct CubicHermiteLocalCoefficients
{

  static constexpr auto makeLocalKeys() {
    if constexpr (dim==1)
      return std::array{
        LocalKey(0, 1, 0),
        LocalKey(0, 1, 1),
        LocalKey(1, 1, 0),
        LocalKey(1, 1, 1)
      };
    if constexpr ((dim==2) and (not reduced))
      return std::array{
        LocalKey(0, 2, 0),
        LocalKey(0, 2, 1),
        LocalKey(0, 2, 2),
        LocalKey(1, 2, 0),
        LocalKey(1, 2, 1),
        LocalKey(1, 2, 2),
        LocalKey(2, 2, 0),
        LocalKey(2, 2, 1),
        LocalKey(2, 2, 2),
        LocalKey(0, 0, 0)
      };
    if constexpr ((dim==2) and (reduced))
      return std::array{
        LocalKey(0, 2, 0),
        LocalKey(0, 2, 1),
        LocalKey(0, 2, 2),
        LocalKey(1, 2, 0),
        LocalKey(1, 2, 1),
        LocalKey(1, 2, 2),
        LocalKey(2, 2, 0),
        LocalKey(2, 2, 1),
        LocalKey(2, 2, 2)
      };
  }

  using LocalKeys = std::decay_t<decltype(makeLocalKeys())>;

public:

  std::size_t size() const
  {
    return localKeys_.size();
  }

  const LocalKey& localKey(std::size_t i) const
  {
    assert( i < localKeys_.size() );
    return localKeys_[i];
  }
private:
  LocalKeys localKeys_ = makeLocalKeys();
};



template<class DF, class RF, unsigned int dim, bool reduced>
class CubicHermiteLocalInterpolation
{
  using ElementJacobianInverse = Dune::FieldMatrix<DF, dim,dim>;
  using LocalSubEntityMeshSize = std::vector<double>;

public:

  template<class Element>
  void bind(const Element& element, const LocalSubEntityMeshSize& localSubEntityMeshSize) {
    localSubEntityMeshSize_ = &localSubEntityMeshSize;
  }

  template<class F, class C>
  void interpolate(const F& f, std::vector<C>& out) const
  {
    using Domain = Dune::FieldVector<DF, dim>;
    auto&& df = derivative(f);
    if constexpr (dim==1)
    {
      out.resize(4);
      out[0] = f(0);
      out[1] = df(0) * (*localSubEntityMeshSize_)[1];
      out[2] = f(1);
      out[3] = df(1) * (*localSubEntityMeshSize_)[3];
    }
    if constexpr (dim==2)
    {
      if constexpr (not reduced)
      {
        out.resize(10);
        out[9] = f(Domain({1.0/3.0,1.0/3.0}));
      }
      if constexpr (reduced)
        out.resize(9);
      auto J0 = df(Domain({0,0}));
      auto J1 = df(Domain({1,0}));
      auto J2 = df(Domain({0,1}));
      out[0] = f(Domain({0,0}));
      out[1] = J0[0][0] * (*localSubEntityMeshSize_)[1];
      out[2] = J0[0][1] * (*localSubEntityMeshSize_)[2];
      out[3] = f(Domain({1,0}));
      out[4] = J1[0][0] * (*localSubEntityMeshSize_)[4];
      out[5] = J1[0][1] * (*localSubEntityMeshSize_)[5];
      out[6] = f(Domain({0,1}));
      out[7] = J2[0][0] * (*localSubEntityMeshSize_)[7];
      out[8] = J2[0][1] * (*localSubEntityMeshSize_)[8];
    }
  }
private:
  const LocalSubEntityMeshSize* localSubEntityMeshSize_;
};



/**
 * \brief Cubic-Hermite shape functions
 *
 * \tparam DF type to represent the field in the domain.
 * \tparam RF type to represent the field in the range.
 * \tparam dim domain dimension
 */
template<class DF, class RF, unsigned int dim, bool reduced>
class CubicHermiteLocalFiniteElement
{
  using LocalBasis = CubicHermiteLocalBasis<DF, RF, dim, reduced>;
  using LocalCoefficients = CubicHermiteLocalCoefficients<dim, reduced>;
  using LocalInterpolation = CubicHermiteLocalInterpolation<DF, RF, dim, reduced>;
  using LocalSubEntityMeshSize = std::vector<double>;

public:

  using Traits = LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>;

  const LocalBasis& localBasis() const {
    return localBasis_;
  }

  const LocalCoefficients& localCoefficients() const {
    return localCoefficients_;
  }

  const LocalInterpolation& localInterpolation() const {
    return localInterpolation_;
  }

  unsigned int size() const {
    return localBasis_.size();
  }

  GeometryType type() const {
    return type_;
  }

  template<class Element, class Mapper, class MeshSizeContainer>
  void bind(const Element& element, const Mapper& mapper, const MeshSizeContainer& subEntityMeshSize) {
    // Update local mesh size cache
    localSubEntityMeshSize_.resize(localCoefficients_.size());
    for(auto i : Dune::range(size()))
    {
      auto localKey = localCoefficients_.localKey(i);
      auto subEntityIndex = mapper.subIndex(element, localKey.subEntity(), localKey.codim());
      localSubEntityMeshSize_[i] = subEntityMeshSize[subEntityIndex];
    }

    localBasis_.bind(element, localSubEntityMeshSize_);
    localInterpolation_.bind(element, localSubEntityMeshSize_);
    type_ = element.type();
  }

private:
  LocalSubEntityMeshSize localSubEntityMeshSize_;
  LocalBasis localBasis_;
  LocalCoefficients localCoefficients_;
  LocalInterpolation localInterpolation_;
  GeometryType type_;
};




} // namespace Impl in Dune::Functions::



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   CubicHermitePreBasis
//   CubicHermiteNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, bool reduced>
class CubicHermiteNode :
  public LeafBasisNode
{
  static const int gridDim = GV::dimension;

  using CubeFiniteElement = Impl::CubicHermiteLocalFiniteElement<typename GV::ctype, double, gridDim, reduced>;
  using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
  using Base = LeafBasisNode;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = CubeFiniteElement;
  using Base::size;

  CubicHermiteNode(const SubEntityMapper& subEntityMapper, const std::vector<double>& subEntityMeshSize) :
    finiteElement_(),
    element_(nullptr),
    subEntityMapper_(&subEntityMapper),
    subEntityMeshSize_(&subEntityMeshSize)
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_.bind(*element_, *subEntityMapper_, *subEntityMeshSize_);
    this->setSize(finiteElement_.size());
  }

protected:

  FiniteElement finiteElement_;
  const Element* element_;
  const SubEntityMapper* subEntityMapper_;
  const std::vector<double>* subEntityMeshSize_;
};



/**
 * \brief Pre-basis for a CubicHermite basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 */
template<typename GV, bool reduced = false>
class CubicHermitePreBasis : public LeafPreBasisMapperMixIn<GV>
{
  using Base = LeafPreBasisMapperMixIn<GV>;
  using Base::mapper_;
  using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

  static constexpr auto cubicHermiteMapperLayout(Dune::GeometryType type, int gridDim) {
    if (type.isVertex())
      return 1 + gridDim;
    if ((type.isTriangle()) and (not reduced))
      return 1;
    else
      return 0;
  }

  static constexpr auto subEntityLayout(Dune::GeometryType type, int gridDim) {
    return (cubicHermiteMapperLayout(type, gridDim) > 0);
  }

public:

  using GridView = typename Base::GridView;
  using Node = CubicHermiteNode<GridView, reduced>;

  CubicHermitePreBasis(const GridView& gv) :
    Base(gv, cubicHermiteMapperLayout),
    subEntityMapper_(gv, subEntityLayout)
  {
    static_assert(GridView::dimension<=2, "CubicHermitePreBasis is only implemented for dim<=2");
    subEntityMeshSize_ = Impl::computeAverageSubEntityMeshSize(subEntityMapper_);
  }

  void update(const GridView& gv)
  {
    Base::update(gv);
    subEntityMapper_.update(gv);
    subEntityMeshSize_ = Impl::computeAverageSubEntityMeshSize(subEntityMapper_);
  }

  Node makeNode() const
  {
    return Node{subEntityMapper_, subEntityMeshSize_};
  }

private:
  std::vector<double> subEntityMeshSize_;
  SubEntityMapper subEntityMapper_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a CubicHermite pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
template<class Dummy=void>
auto cubicHermite()
{
  return [](const auto& gridView) {
    return CubicHermitePreBasis<std::decay_t<decltype(gridView)>>(gridView);
  };
}

/**
 * \brief Create a pre-basis factory that can create a CubicHermite pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
template<class Dummy=void>
auto reducedCubicHermite()
{
  return [](const auto& gridView) {
    return CubicHermitePreBasis<std::decay_t<decltype(gridView)>, true>(gridView);
  };
}

} // end namespace BasisFactory




/** \brief CubicHermite basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using CubicHermiteBasis = DefaultGlobalBasis<CubicHermitePreBasis<GV>>;



/** \brief CubicHermite basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using ReducedCubicHermiteBasis = DefaultGlobalBasis<CubicHermitePreBasis<GV, true>>;



} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CUBICHERMITEBASIS_HH
