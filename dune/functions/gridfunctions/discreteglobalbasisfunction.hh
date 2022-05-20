// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH

#include <memory>

#include <dune/common/typetraits.hh>

#include <dune/typetree/treecontainer.hh>

#include <dune/functions/functionspacebases/hierarchicnodetorangemap.hh>
#include <dune/functions/functionspacebases/flatvectorview.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/backends/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

namespace Dune {
namespace Functions {



/**
 * \brief A grid function induced by a global basis and a coefficient vector.
 *
 * \ingroup FunctionImplementations
 *
 * This implements the grid function interface by combining a given global
 * basis and a coefficient vector.
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
 * \tparam B Type of global basis
 * \tparam V Type of coefficient vectors
 * \tparam NTRE Type of node-to-range-entry-map that associates each leaf node in the local ansatz subtree with an entry in the range type
 * \tparam R Range type of this function
 */
template<typename B, typename V,
  typename NTRE = HierarchicNodeToRangeMap,
  typename R = typename V::value_type>
class DiscreteGlobalBasisFunction
{
public:
  using Basis = B;
  using Vector = V;

  // In order to make the cache work for proxy-references
  // we have to use AutonomousValue<T> instead of std::decay_t<T>
  using Coefficient = Dune::AutonomousValue<decltype(std::declval<Vector>()[std::declval<typename Basis::MultiIndex>()])>;

  using GridView = typename Basis::GridView;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Tree = typename Basis::LocalView::Tree;
  using NodeToRangeEntry = NTRE;

  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = R;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  using Traits = Imp::GridFunctionTraits<Range(Domain), EntitySet, DefaultDerivativeTraits, 16>;

  class LocalFunction
  {
    using LocalView = typename Basis::LocalView;
    using size_type = typename Tree::size_type;

    template<class Node>
    using LocalBasisRange = typename Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;

    template<class Node>
    using NodeData = typename std::vector<LocalBasisRange<Node>>;

    using PerNodeEvaluationBuffer = typename TypeTree::TreeContainer<NodeData,Tree>;

  public:

    using GlobalFunction = DiscreteGlobalBasisFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using Element = GlobalFunction::Element;

    //! Create a local-function from the associated grid-function
    LocalFunction(const DiscreteGlobalBasisFunction& globalFunction)
      : globalFunction_(&globalFunction)
      , localView_(globalFunction.basis().localView())
      , evaluationBuffer_(localView_.tree())
    {
      localDoFs_.reserve(localView_.maxSize());
    }

    /**
     * \brief Copy-construct the local-function.
     *
     * This copy-constructor copies the cached local DOFs only
     * if the `other` local-function is bound to an element.
     **/
    LocalFunction(const LocalFunction& other)
      : globalFunction_(other.globalFunction_)
      , localView_(other.localView_)
      , evaluationBuffer_(localView_.tree())
    {
      localDoFs_.reserve(localView_.maxSize());
      if (bound())
        localDoFs_ = other.localDoFs_;
    }

    /**
     * \brief Copy-assignment of the local-function.
     *
     * Assign all members from `other` to `this`, except the
     * local DOFs. Those are copied only if the `other`
     * local-function is bound to an element.
     **/
    LocalFunction& operator=(const LocalFunction& other)
    {
      globalFunction_ = other.globalFunction_;
      localView_ = other.localView_;
      if (bound())
        localDoFs_ = other.localDoFs_;
      return *this;
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before `operator()`
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element)
    {
      localView_.bind(element);
      // Use cache of full local view size. For a subspace basis,
      // this may be larger than the number of local DOFs in the
      // tree. In this case only cache entries associated to local
      // DOFs in the subspace are filled. Cache entries associated
      // to local DOFs which are not contained in the subspace will
      // not be touched.
      //
      // Alternatively one could use a cache that exactly fits
      // the size of the tree. However, this would require to
      // subtract an offset from localIndex(i) on each cache
      // access in operator().
      localDoFs_.resize(localView_.size());
      for (size_type i = 0; i < localView_.tree().size(); ++i)
      {
        // For a subspace basis the index-within-tree i
        // is not the same as the localIndex withn the
        // full local view.
        size_t localIndex = localView_.tree().localIndex(i);
        localDoFs_[localIndex] = globalFunction_->dofs()[localView_.index(localIndex)];
      }
    }

    //! Unbind the local-function.
    void unbind()
    {
      localView_.unbind();
    }

    //! Check if LocalFunction is already bound to an element.
    bool bound() const
    {
      return localView_.bound();
    }

    /**
     * \brief Evaluate this local-function in coordinates `x` in the bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     */
    Range operator()(const Domain& x) const
    {
      Range y;
      istlVectorBackend(y) = 0;

      TypeTree::forEachLeafNode(localView_.tree(), [&](auto&& node, auto&& treePath) {
        const auto& nodeToRangeEntry = globalFunction_->nodeToRangeEntry();
        const auto& fe = node.finiteElement();
        const auto& localBasis = fe.localBasis();
        auto& shapeFunctionValues = evaluationBuffer_[treePath];

        localBasis.evaluateFunction(x, shapeFunctionValues);

        // Get range entry associated to this node
        auto re = flatVectorView(nodeToRangeEntry(node, treePath, y));

        for (size_type i = 0; i < localBasis.size(); ++i)
        {
          // Get coefficient associated to i-th shape function
          auto c = flatVectorView(localDoFs_[node.localIndex(i)]);

          // Get value of i-th shape function
          auto v = flatVectorView(shapeFunctionValues[i]);

          // Notice that the range entry re, the coefficient c, and the shape functions
          // value v may all be scalar, vector, matrix, or general container valued.
          // The matching of their entries is done via the multistage procedure described
          // in the class documentation of DiscreteGlobalBasisFunction.
          auto&& dimC = c.size();
          auto dimV = v.size();
          assert(dimC*dimV == re.size());
          for(size_type j=0; j<dimC; ++j)
          {
            auto&& c_j = c[j];
            for(size_type k=0; k<dimV; ++k)
              re[j*dimV + k] += c_j*v[k];
          }
        }
      });

      return y;
    }

    //! Return the element the local-function is bound to.
    const Element& localContext() const
    {
      return localView_.element();
    }

    //! Not implemented.
    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:

    const DiscreteGlobalBasisFunction* globalFunction_;
    LocalView localView_;
    mutable PerNodeEvaluationBuffer evaluationBuffer_;
    std::vector<Coefficient> localDoFs_;
  };

  //! Create a grid-function, by wrapping the arguments in `std::shared_ptr`.
  template<class B_T, class V_T, class NTRE_T>
  DiscreteGlobalBasisFunction(B_T && basis, V_T && coefficients, NTRE_T&& nodeToRangeEntry) :
    entitySet_(basis.gridView()),
    basis_(wrap_or_move(std::forward<B_T>(basis))),
    coefficients_(wrap_or_move(std::forward<V_T>(coefficients))),
    nodeToRangeEntry_(wrap_or_move(std::forward<NTRE_T>(nodeToRangeEntry)))
  {}

  //! Create a grid-function, by moving the arguments in `std::shared_ptr`.
  DiscreteGlobalBasisFunction(std::shared_ptr<const Basis> basis, std::shared_ptr<const V> coefficients, std::shared_ptr<const NodeToRangeEntry> nodeToRangeEntry) :
    entitySet_(basis->gridView()),
    basis_(basis),
    coefficients_(coefficients),
    nodeToRangeEntry_(nodeToRangeEntry)
  {}

  //! Return a const reference to the stored basis.
  const Basis& basis() const
  {
    return *basis_;
  }

  //! Return the coefficients of this discrete function by reference.
  const V& dofs() const
  {
    return *coefficients_;
  }

  //! Return the stored node-to-range map.
  const NodeToRangeEntry& nodeToRangeEntry() const
  {
    return *nodeToRangeEntry_;
  }

  //! Not implemented.
  Range operator() (const Domain& x) const
  {
    // TODO: Implement this using hierarchic search
    DUNE_THROW(NotImplemented,"not implemented");
  }

  //! Not implemented.
  friend typename Traits::DerivativeInterface derivative(const DiscreteGlobalBasisFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  /**
   * \brief Construct local function from a DiscreteGlobalBasisFunction.
   *
   * The obtained a local-function the satisfies the concept
   * `Dune::Functions::Concept::LocalFunction`. It must be bound
   * to an entity from the entity set of the DiscreteGlobalBasisFunction
   * before it can be used.
   *
   * Notice that the local-function stores a reference to the
   * global DiscreteGlobalBasisFunction. Hence calling any method
   * of the local-function after the DiscreteGlobalBasisFunction
   * exceeded its life time leads to undefined behavior.
   */
  friend LocalFunction localFunction(const DiscreteGlobalBasisFunction& t)
  {
    return LocalFunction(t);
  }

  //! Get associated set of entities the local-function can be bound to.
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:

  EntitySet entitySet_;
  std::shared_ptr<const Basis> basis_;
  std::shared_ptr<const V> coefficients_;
  std::shared_ptr<const NodeToRangeEntry> nodeToRangeEntry_;
};



/**
 * \brief Construction of local-functions from a temporary DiscreteGlobalBasisFunction (forbidden).
 *
 * Since a DiscreteGlobalBasisFunction::LocalFunction stores a reference
 * to the global DiscreteGlobalBasisFunction its life time is bound to
 * the latter. Hence construction from a temporary DiscreteGlobalBasisFunction
 * would lead to a dangling reference and is thus forbidden/deleted.
 */
template<typename... TT>
void localFunction(DiscreteGlobalBasisFunction<TT...>&& t) = delete;


/**
 * \brief Generate a DiscreteGlobalBasisFunction.
 *
 * \ingroup FunctionImplementations
 *
 * Create a new DiscreteGlobalBasisFunction by wrapping the vector in a
 * VectorBackend that allows the hierarchic resize and multi-index access in
 * the DiscreteGlobalBasisFunction, if the vector does not yet fulfill the
 * \ref ConstVectorBackend concept.
 *
 * \tparam R  The range type this grid-function should represent when seen as
 *            a mapping `R(Domain)` with `Domain` the global coordinates of the
 *            associated GridView. This must be compatible with the basis and
 *            coefficients. See the documentation of \ref DiscreteGlobalBasisFunction
 *            for more details.
 *
 * \param basis  The global basis or subspace basis associated with this
 *               grid-function
 * \param vector The coefficient vector to use in combination with the `basis`.
 *
 * \relatesalso DiscreteGlobalBasisFunction
 **/
template<typename R, typename B, typename V>
auto makeDiscreteGlobalBasisFunction(B&& basis, V&& vector)
{
  using Basis = std::decay_t<B>;
  using NTREM = HierarchicNodeToRangeMap;

  // Small helper functions to wrap vectors using istlVectorBackend
  // if they do not already satisfy the VectorBackend interface.
  auto toConstVectorBackend = [&](auto&& v) -> decltype(auto) {
    if constexpr (models<Concept::ConstVectorBackend<Basis>, decltype(v)>()) {
      return std::forward<decltype(v)>(v);
    } else {
      return istlVectorBackend(v);
    }
  };

  using Vector = std::decay_t<decltype(toConstVectorBackend(std::forward<V>(vector)))>;
  return DiscreteGlobalBasisFunction<Basis, Vector, NTREM, R>(
      std::forward<B>(basis),
      toConstVectorBackend(std::forward<V>(vector)),
      HierarchicNodeToRangeMap());
}



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
