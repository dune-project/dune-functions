// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH

#include <array>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/common/hybridutilities.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/std/optional.hh>
#include <dune/common/typeutilities.hh>
#include <dune/typetree/childextraction.hh>

#include <dune/functions/backends/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/treedata.hh>
#include <dune/functions/common/wrapownshare.hh>
#include <dune/functions/functionspacebases/flatvectorview.hh>
#include <dune/functions/functionspacebases/hierarchicnodetorangemap.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>


namespace Dune {
namespace Functions {



/**
 * \brief A grid function induced by a global basis and a coefficient vector
 *
 * \ingroup FunctionImplementations
 *
 * This implements the grid function interface by combining a given global
 * basis and a coefficient vector. The part of the spanned space that should
 * be covered by the function is determined by a tree path that specifies the
 * corresponding local ansatz tree.
 *
 * This class supports mapping of subtrees to multi-component ranges,
 * vector-valued shape functions, and implicit product spaces given
 * by vector-valued coefficients. The mapping of these to the range
 * type is done via the following multistage procedure:
 *
 * 1.Each leaf node N in the local ansatz subtree is associated to an entry
 *   RE of the range-type via the given node-to-range-entry-map.
 *
 * Now let C be the coefficient block for a single basis function and
 * V the value of this basis function at the evaluation point. Notice
 * that both may be scalar, vector, matrix, or general container valued.
 *
 * 2.Each entry of C is associated with a flat index j via flatVectorView.
 *   This is normally a lexicographic index. The total scalar dimension according
 *   to those flat indices is dim(C).
 * 3.Each entry of V is associated with a flat index k via flatVectorView.
 *   This is normally a lexicographic index. The total scalar dimension according
 *   to those flat indices dim(V).
 * 4.Each entry of RE is associated with a flat index k via flatVectorView.
 *   This is normally a lexicographic index. The total scalar dimension according
 *   to those flat indices dim(RE).
 * 5.Via those flat indices we now interpret C,V, and RE as vectors and compute the diadic
 *   product (C x V). The entries of this product are mapped to the flat indices for
 *   RE lexicographically. I.e. we set
 *
 *     RE[j*dim(V)+k] = C[j] * V[k]
 *
 * Hence the range entry RE must have dim(RE) = dim(C)*dim(V).
 *
 * \tparam B Type of lobal basis
 * \tparam TP Type of tree path specifying the requested subtree of ansatz functions
 * \tparam V Type of coefficient vectors
 * \tparam NTRE Type of node-to-range-entry-map that associates each leaf node in the local ansatz subtree with an entry in the range type
 * \tparam R Range type of this function
 */
template<typename B, typename TP, typename V,
  typename NTRE = HierarchicNodeToRangeMap,
  typename R = typename V::value_type>
class DiscreteGlobalBasisFunction
{
public:
  using Basis = B;
  using TreePath = TP;
  using Vector = V;

  using GridView = typename Basis::GridView;
  using LocalView = typename Basis::LocalView;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Tree = typename LocalView::Tree;
  using SubTree = typename TypeTree::ChildForTreePath<Tree, TreePath>;
  using NodeToRangeEntry = NTRE;

  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = R;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  // some internal template defaults
  template <class Signature>
  using DerivativeTraits = DefaultDerivativeTraits<Signature>;
  static constexpr size_t bufferSize = 16;
  template <class Signature>
  using GridFunctionTraits = Imp::GridFunctionTraits<Signature, EntitySet, DerivativeTraits, bufferSize>;

  // function signature and traits
  using Signature = Range(Domain);
  using Traits = GridFunctionTraits<Signature>;

  // derivative signature and traits
  using DSignature = typename SignatureTraits<Signature>::template DerivativeSignature<DerivativeTraits>;
  using DTraits = GridFunctionTraits<DSignature>;

  // internal node traits
  template <class Node>
  struct NodeTraits {
    using FETraits = typename Node::FiniteElement::Traits;
    using LBTraits = typename FETraits::LocalBasisType::Traits;
    static constexpr size_t dimDomain = LBTraits::dimDomain;
    static constexpr size_t dimRange = LBTraits::dimRange;
    using Range = typename LBTraits::RangeType;
    using Jacobian = typename LBTraits::JacobianType;
  };

  class LocalFunction
  {
    template<class Node>
    using NodeData = typename std::vector<typename NodeTraits<Node>::Range>;

    using ShapeFunctionValueContainer = TreeData<SubTree, NodeData, true>;

    struct LocalEvaluateVisitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {
      LocalEvaluateVisitor(const LocalDomain& x, Range& y, const LocalView& localView, const Vector& coefficients, const NodeToRangeEntry& nodeToRangeEntry, ShapeFunctionValueContainer& shapeFunctionValueContainer):
        x_(x),
        y_(y),
        localView_(localView),
        coefficients_(coefficients),
        nodeToRangeEntry_(nodeToRangeEntry),
        shapeFunctionValueContainer_(shapeFunctionValueContainer)
      {}

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath treePath)
      {
        auto&& fe = node.finiteElement();
        auto&& localBasis = fe.localBasis();

        auto& shapeFunctionValues = shapeFunctionValueContainer_[node];
        localBasis.evaluateFunction(x_, shapeFunctionValues);

        // Get range entry associated to this node
        auto re = flatVectorView(nodeToRangeEntry_(node, treePath, y_));

        for (size_t i = 0; i < localBasis.size(); ++i)
        {
          auto&& multiIndex = localView_.index(node.localIndex(i));

          // Get coefficient associated to i-th shape function
          auto c = flatVectorView(coefficients_[multiIndex]);

          // Get value of i-th shape function
          auto v = flatVectorView(shapeFunctionValues[i]);

          // Notice that the range entry re, the coefficient c, and the shape functions
          // value v may all be scalar, vector, matrix, or general container valued.
          // The matching of their entries is done via the multistage procedure described
          // in the class documentation of DiscreteGlobalBasisFunction.
          auto&& dimC = c.size();
          auto dimV = v.size();
          assert(dimC*dimV == re.size());
          for(size_t j=0; j<dimC; ++j)
          {
            auto&& c_j = c[j];
            for(size_t k=0; k<dimV; ++k)
              re[j*dimV + k] += c_j*v[k];
          }
        }
      }

    private:
      const LocalDomain& x_;
      Range& y_;
      const LocalView& localView_;
      const Vector& coefficients_;
      const NodeToRangeEntry& nodeToRangeEntry_;
      ShapeFunctionValueContainer& shapeFunctionValueContainer_;
    };


  public:
    using GlobalFunction = DiscreteGlobalBasisFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using Element = GlobalFunction::Element;

    LocalFunction(const DiscreteGlobalBasisFunction& globalFunction)
      : globalFunction_(globalFunction)
      , localView_(globalFunction.basis().localView())
    {
      init();
    }

    LocalFunction(const LocalFunction& other)
      : globalFunction_(other.globalFunction_)
      , localView_(other.localView_)
    {
      init();
    }

    LocalFunction operator=(const LocalFunction& other)
    {
      globalFunction_ = other.globalFunction_;
      localView_ = other.localView_;
      init();
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluate()
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element) { localView_.bind(element); }

    void unbind() { localView_.unbind(); }

    /**
     * \brief Check if LocalFunction is already bound to an element.
     */
    bool bound() const { return localView_.isBound(); }

    /**
     * \brief Return element that LocalFunction is bound to.
     */
    const Element& localContext() const { return localView_.element(); }

    /**
     * \brief Evaluate LocalFunction at bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     */
    Range operator()(const Domain& x) const
    {
      auto y = Range(0);

      LocalEvaluateVisitor localEvaluateVisitor(x, y, localView_, globalFunction_.dofs(), globalFunction_.nodeToRangeEntry(), shapeFunctionValueContainer_);
      TypeTree::applyToTree(*subTree_, localEvaluateVisitor);

      return y;
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction&)
    {
      DUNE_THROW(NotImplemented, "Derivative of LocalFunction of DiscreteGlobalBasisFunction.");
    }

  protected:
    const DiscreteGlobalBasisFunction& globalFunction_;
    LocalView localView_;

    mutable ShapeFunctionValueContainer shapeFunctionValueContainer_;
    const SubTree* subTree_;

  private:
    void init() {
      // Here we assume that the tree can be accessed, traversed,
      // and queried for tree indices even in unbound state.
      using namespace TypeTree;
      subTree_ = &child(localView_.tree(), globalFunction_.treePath());
      shapeFunctionValueContainer_.init(*subTree_);
    }
  };

  template<class B_T, class V_T, class NTRE_T>
  DiscreteGlobalBasisFunction(B_T && basis, const TreePath& treePath, V_T && coefficients, NTRE_T&& nodeToRangeEntry) :
    entitySet_(basis.gridView()),
    basis_(wrap_own_share<const B>(std::forward<B_T>(basis))),
    treePath_(treePath),
    coefficients_(wrap_own_share<const V>(std::forward<V_T>(coefficients))),
    nodeToRangeEntry_(wrap_own_share<const NTRE>(std::forward<NTRE_T>(nodeToRangeEntry)))
  {}

  // getters
  const auto& basis() const { return *basis_; }
  const auto& treePath() const { return treePath_; }
  const auto& dofs() const { return *coefficients_; }
  const auto& nodeToRangeEntry() const { return *nodeToRangeEntry_; }
  const auto& entitySet() const { return entitySet_; }

  // TODO: Implement this using hierarchic search
  Range operator() (const Domain&) const
  {
    DUNE_THROW(NotImplemented, "Evaluation of DiscreteGlobalBasisFunction.");
  }

  /**
   * @brief Derivative of DiscreteGlobalFunction
   */
  struct Derivative {

    Derivative(const DiscreteGlobalBasisFunction& f)
        : f_(f) {}

    const EntitySet& entitySet() const { return f_.entitySet(); }

    typename DTraits::Range operator()(const Domain&) const {
      DUNE_THROW(NotImplemented,
                 "Evaluation of Derivative of DiscreteGlobalBasisFunction.");
    }

    /**
     * @brief LocalFunction of Derivative of DiscreteGlobalBasisFunction
     */
    struct LocalDerivative {

      using Geometry = std::decay_t<typename Element::Geometry>;
      using JInverseT = typename Geometry::JacobianInverseTransposed;

      template <class Node>
      using NodeData = std::vector<typename NodeTraits<Node>::Jacobian>;
      using ShapeFunctionJacobianContainer = TreeData<SubTree, NodeData, true>;

      using LocalDRange = typename DTraits::LocalFunctionTraits::Range;

      struct LocalDerivativeVisitor
          : TypeTree::TreeVisitor
          , TypeTree::StaticTraversal {

        // Note: we assume that the NodeToRangeEntry can also be used as
        // a NodeToJacobianEntry. See all lines tagged with (NTJE)
        using NodeToJacobianEntry = NodeToRangeEntry;

        LocalDerivativeVisitor(
            const LocalDomain& x, LocalDRange& y, const LocalView& localView,
            const Vector& coefficients,
            const NodeToJacobianEntry& nodeToJacobianEntry, // see (NTJE)
            ShapeFunctionJacobianContainer& shapeFunctionJacobianContainer,
            JInverseT jacobianInverseTransposed)
            : x_(x)
            , y_(y)
            , localView_(localView)
            , coefficients_(coefficients)
            , nodeToJacobianEntry_(nodeToJacobianEntry) // see (NTJE)
            , shapeFunctionJacobianContainer_(shapeFunctionJacobianContainer)
            , jInverseT_(std::move(jacobianInverseTransposed)) {}

        template <typename Node, typename TreePath>
        void leaf(const Node& node, TreePath treePath) {

          auto&& fe = node.finiteElement();
          auto&& localBasis = fe.localBasis();

          // TODO we can cache information here for power-nodes.
          auto&& shapeFunctionJacobians = shapeFunctionJacobianContainer_[node];
          localBasis.evaluateJacobian(x_, shapeFunctionJacobians);

          /// The following approach first computes the
          /// (1)   referenceJacobian = sum_i c_i * jacobian(localBasis_i)
          /// and subsequently transforms them
          /// (2)        y = referenceJacobian * jacobianInverseTransposed

          // Notice that the range entry re, the coefficient c, and the shape
          // functions value v may all be scalar, vector, matrix, or general
          // container valued. The matching of their entries is done via the
          // multistage procedure described in the class documentation of
          // DiscreteGlobalBasisFunction.

          /// part (1)
          // Since the coefficients may be non-scalar, we need a way to
          // store the respective results. We use an array to store
          // the referenceJacobian_k for each (c_i)_k where j is a flat index
          // for whatever (statically sized) type (c_i) is.
          constexpr size_t dimC = 1;
          // TODO generalize to non-scalar coefficients. A first improvement
          // might deduce a (uniform) size. A second improvement might deduce
          // node-dependet sizes, e.g. using a TreeData mechanism.
          using Jacobian = typename NodeTraits<Node>::Jacobian;
          using ReferenceJ = std::array<Jacobian, dimC>;
          ReferenceJ referenceJData;
          for (size_t k = 0; k < dimC; ++k)
            referenceJData[k] = 0;
          auto referenceJ = flatVectorView(referenceJData);
          size_t dimReferenceJ = referenceJ.size();

          // Compute the reference jacobians
          for (size_t i = 0; i < localBasis.size(); ++i) {
            auto&& multiIndex = localView_.index(node.localIndex(i));
            // get i-th coefficient and jacobian of i-th local basis
            auto&& c_i = flatVectorView(coefficients_[multiIndex]);
            auto&& j_i = flatVectorView(shapeFunctionJacobians[i]);
            assert(dimC == c_i.size() && "Implement non-scalar coefficients.");
            size_t dimJ = j_i.size();
            assert(dimC * dimJ == dimReferenceJ);
            for (size_t k = 0; k < dimC; ++k)
              for (size_t d = 0; d < dimJ; ++d)
                referenceJ[k * dimJ + d] += c_i[k] * j_i[d];
          }

          /// part (2)
          // Provide flat access
          auto&& re = nodeToJacobianEntry_(node, treePath, y_); // see (NTJE)
          auto&& y = flatVectorView(re);
          size_t dimRange = y.size();
          auto&& jInverseT = flatVectorView(jInverseT_);

          // Assert some size connections
          size_t dimLocalReferenceBasisDomain = NodeTraits<Node>::dimDomain;
          size_t dimLocalReferenceBasisRange = NodeTraits<Node>::dimRange;
          using Geometry = std::decay_t<decltype(node.element().geometry())>;
          assert(Geometry::mydimension == dimLocalReferenceBasisDomain);
          size_t dimJIT = jInverseT.size();
          assert(dimJIT == Geometry::mydimension * Geometry::coorddimension);
          size_t rowsRange = dimC * dimLocalReferenceBasisRange;
          size_t colsRange = Geometry::coorddimension;
          assert(dimRange == rowsRange * colsRange);

          /*
           * y = referenceGradients * jacobianInverse
           *
           *             and in terms of dimensions
           * [rowsRange x colsRange] =
           *         [dimC*dimLocalBasisRange x dimLocalReferenceBasisDomain]
           *  times  [dimLocalReferenceBasisDomain x Geometry::coordimension]
           *
           */
          // TODO check ordering when non-scalar coefficients are used.
          for (size_t row = 0; row < rowsRange; ++row) {
            for (size_t col = 0; col < colsRange; ++col) {
              auto&& range = y[row * colsRange + col];
              for (size_t l = 0; l < dimLocalReferenceBasisDomain; ++l) {
                auto&& ref = referenceJ[row * dimLocalReferenceBasisDomain + l];
                // wiggle indices to untranspose jacobianInverseTransposed
                auto&& jit = jInverseT[col * dimLocalReferenceBasisDomain + l];
                range += ref * jit;
              }
            }
          }
        }

      private:
        const LocalDomain& x_;
        LocalDRange& y_;
        const LocalView& localView_;
        const Vector& coefficients_;
        const NodeToRangeEntry& nodeToJacobianEntry_; // see (NTJE)
        ShapeFunctionJacobianContainer& shapeFunctionJacobianContainer_;
        JInverseT jInverseT_;
      };

      template <class D_, disableCopyMove<LocalDerivative, D_> = 0>
      LocalDerivative(D_&& derivative)
          : derivative_(wrap_own_share<const Derivative>(std::forward<D_>(derivative)))
          , localView_(derivative_->basis().localView()) {
        init();
      }

      LocalDerivative(const LocalDerivative& other)
          : derivative_(other.derivative_)
          , localView_(other.localView_) {
        init();
      }

      LocalDerivative operator=(const LocalDerivative& other) {
        derivative_ = other.dgbfDerivative;
        localView_ = other.localView_;
        init();
      }

      void bind(const Element& e) {
        localView_.bind(e);
        geometry_.emplace(e.geometry());
      }
      void unbind() {
        geometry_.reset();
        localView_.unbind();
      }

      const Element& localContext() const { return localView_.element(); }

      LocalDRange operator()(const Domain& x) const {
        LocalDRange y;
        y = 0;

        LocalDerivativeVisitor localDerivativeVisitor(
            x, y, localView_, derivative_->dofs(),
            derivative_->nodeToRangeEntry(),
            shapeFunctionJacobianContainer_,
            geometry_.value().jacobianInverseTransposed(x));
        TypeTree::applyToTree(*subTree_, localDerivativeVisitor);

        return y;
      }

      friend typename DTraits::LocalFunctionTraits::DerivativeInterface
      derivative(const LocalDerivative&) {
        DUNE_THROW(NotImplemented, "Derivative of LocalFunction of Derivative "
                                   "of DiscreteGlobalBasisFunction.");
      }

    private:
      std::shared_ptr<const Derivative> derivative_;
      LocalView localView_;
      Std::optional<Geometry> geometry_;
      mutable ShapeFunctionJacobianContainer shapeFunctionJacobianContainer_;
      const SubTree* subTree_;

      void init() {
        // Here we assume that the tree can be accessed, traversed,
        // and queried for tree indices even in unbound state.
        using namespace TypeTree;
        subTree_ = &child(localView_.tree(), derivative_->treePath());
        shapeFunctionJacobianContainer_.init(*subTree_);
      }
    };

    /**
     * \brief Construct local function from a Derivative of a DiscreteGlobalBasisFunction
     *
     * \ingroup FunctionImplementations
     *
     * The obtained local function satisfies the concept
     * `Dune::Functions::Concept::LocalFunction` and must be bound
     * to an entity from the entity set of the DiscreteGlobalBasisFunction (or,
     * equivalently, its Derivative) before it can be used.
     *
     * Notice that the local function shares ownership for the Derivative of the
     * global DiscreteGlobalBasisFunction, i.e. life time of the Derivative is
     * extended appropriately (if constructed correctly).
     * The underlying DiscreteGlobalBasisFunction needs to live though, since
     * its derivative stores a reference (currently). Hence calling any method
     * of this local function after the DiscreteGlobalBasisFunction exceeded its
     * life time leads to undefined behavior.
     */
    friend auto localFunction(const Derivative& d) { return LocalDerivative(d); }
    friend auto localFunction(Derivative&& d) { return LocalDerivative(std::move(d)); }

    friend typename DTraits::DerivativeInterface
    derivative(const Derivative&) {
      DUNE_THROW(NotImplemented,
                 "Derivative of Derivative of DiscreteGlobalBasisFunction.");
    }

    // getters
    const auto& basis() const { return f_.basis(); }
    const auto& treePath() const { return f_.treePath(); }
    const auto& dofs() const { return f_.dofs(); }
    const auto& nodeToRangeEntry() const { return f_.nodeToRangeEntry(); }

  private:
    const DiscreteGlobalBasisFunction& f_;
  };

  /**
   * \brief Construct Derivative of a DiscreteGlobalBasisFunction
   *
   * \ingroup FunctionImplementations
   *
   * The obtained global function can be used to create local functions
   * that yield local derivative evaluations.
   *
   */
  friend auto derivative(const DiscreteGlobalBasisFunction& t) {
    return Derivative(t);
  }

  /**
   * \brief Construct local function from a DiscreteGlobalBasisFunction
   *
   * \ingroup FunctionImplementations
   *
   * The obtained local function satisfies the concept
   * `Dune::Functions::Concept::LocalFunction` and must be bound
   * to an entity from the entity set of the DiscreteGlobalBasisFunction
   * before it can be used.
   *
   * Notice that the local function stores a reference to the
   * global DiscreteGlobalBasisFunction. Hence calling any method
   * of the local function after the DiscreteGlobalBasisFunction
   * exceeded its life time leads to undefined behavior.
   */
  friend LocalFunction localFunction(const DiscreteGlobalBasisFunction& t)
  {
    return LocalFunction(t);
  }

private:
  EntitySet entitySet_;
  std::shared_ptr<const B> basis_;
  const TreePath treePath_;
  std::shared_ptr<const V> coefficients_;
  std::shared_ptr<const NTRE> nodeToRangeEntry_;
};



/**
 * \brief Construction of local functions from a temporary DiscreteGlobalBasisFunction (forbidden)
 *
 * Since a DiscreteGlobalBasisFunction::LocalFunction stores a reference
 * to the global DiscreteGlobalBasisFunction its life time is bound to
 * the latter. Hence construction from a temporary DiscreteGlobalBasisFunction
 * would lead to a dangling reference and is thus forbidden/deleted.
 *
 * \ingroup FunctionImplementations
 */
template<typename... TT>
void localFunction(DiscreteGlobalBasisFunction<TT...>&& t) = delete;



template<typename R, typename B, typename TP, typename V>
auto makeDiscreteGlobalBasisFunction(B&& basis, const TP& treePath, V&& vector)
{
  using Basis = std::decay_t<B>;
  using NTREM = HierarchicNodeToRangeMap;

  // Small helper functions to wrap vectors using istlVectorBackend
  // if they do not already satisfy the VectorBackend interface.
  auto toConstVectorBackend = [&](auto&& v) -> decltype(auto) {
    return Dune::Hybrid::ifElse(models<Concept::ConstVectorBackend<Basis>, decltype(v)>(),
    [&](auto id) -> decltype(auto) {
      return std::forward<decltype(v)>(v);
    }, [&](auto id) -> decltype(auto) {
      return istlVectorBackend(v);
    });
  };

  using Vector = std::decay_t<decltype(toConstVectorBackend(std::forward<V>(vector)))>;
  return DiscreteGlobalBasisFunction<Basis, TP, Vector, NTREM, R>(
      std::forward<B>(basis),
      treePath,
      toConstVectorBackend(std::forward<V>(vector)),
      HierarchicNodeToRangeMap());
}



template<typename R, typename B, typename V>
auto makeDiscreteGlobalBasisFunction(B&& basis, V&& vector)
{
  return makeDiscreteGlobalBasisFunction<R>(std::forward<B>(basis), TypeTree::hybridTreePath(), std::forward<V>(vector));
}



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
