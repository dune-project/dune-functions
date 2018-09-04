// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH

#include <memory>

#include <dune/common/shared_ptr.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/functions/functionspacebases/hierarchicnodetorangemap.hh>
#include <dune/functions/functionspacebases/flatvectorview.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/common/treedata.hh>
#include <dune/functions/common/wrapownshare.hh>
#include <dune/functions/backends/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

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

  // internal node traits
  template <class Node>
  struct NodeTraits {
    using FETraits = typename Node::FiniteElement::Traits;
    using LBTraits = typename FETraits::LocalBasisType::Traits;
    using Range = typename LBTraits::RangeType;
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

        auto&& shapeFunctionValues = shapeFunctionValueContainer_[node];
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

  friend typename Traits::DerivativeInterface derivative(const DiscreteGlobalBasisFunction&)
  {
    DUNE_THROW(NotImplemented, "Derivative of DiscreteGlobalBasisFunction.");
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
