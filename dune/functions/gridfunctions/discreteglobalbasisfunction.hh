// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH

#include <memory>

#include <dune/common/shared_ptr.hh>

#include <dune/localfunctions/finiteelementtypes/functionspacetypes.hh>

#include <dune/functions/functionspacebases/defaultnodetorangemap.hh>
#include <dune/functions/functionspacebases/flatvectorbackend.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/common/treedata.hh>



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
 * 2.Each entry of C is associated with a flat index j via FlatVectorBackend.
 *   This is normally a lexicographic index. The total scalar dimension according
 *   to those flat indices is dim(C).
 * 3.Each entry of V is associated with a flat index k via FlatVectorBackend.
 *   This is normally a lexicographic index. The total scalar dimension according
 *   to those flat indices dim(V).
 * 4.Each entry of RE is associated with a flat index k via FlatVectorBackend.
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
  typename NTRE = DefaultNodeToRangeMap<typename std::decay<decltype(std::declval<B>().localView().tree().child(std::declval<TP>()))>::type>,
  typename R = typename V::value_type>
class DiscreteGlobalBasisFunction
{
public:
  using Basis = B;
  using TreePath = TP;
  using Vector = V;

  using GridView = typename Basis::GridView;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Tree = typename Basis::LocalView::Tree;
  using SubTree = typename TypeTree::ChildForTreePath<Tree, TreePath>;
  using NodeToRangeEntry = NTRE;

  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = R;
  // How it should be (not working for Range = double):
  //using RangeField = typename Range::value_type;
  //using JacobianRange = FieldVector< FieldVector<RangeField, Domain::dimension>, Range::dimension>;
  using RangeField = double;
  // For vectorial ranges, this works only if Range::dimension == Domain::dimension, otherwise the assert
  // in the evaluation method for the derivative throws an error.
  using JacobianRange = FieldVector< Range, Domain::dimension>;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  using Traits = Imp::GridFunctionTraits<Range(Domain), EntitySet, DefaultDerivativeTraits, 16>;

  class LocalFunction
  {
    using LocalBasisView = typename Basis::LocalView;
    using LocalIndexSet = typename Basis::LocalIndexSet;
    using size_type = typename SubTree::size_type;

    template<class Node>
    using LocalBasisRange = typename Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
    template<class Node>
    using LocalBasisJacobian = typename Node::FiniteElement::Traits::LocalBasisType::Traits::JacobianType;

    template<class Node>
    using NodeData = typename std::vector<LocalBasisRange<Node>>;
    template<class Node>
    using NodeDataJacobian = typename std::vector<LocalBasisJacobian<Node>>;

    using ShapeFunctionValueContainer = TreeData<SubTree, NodeData, true>;
    using ShapeFunctionJacobianContainer = TreeData<SubTree, NodeDataJacobian, true>;

    struct LocalEvaluateVisitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {
      LocalEvaluateVisitor(const LocalDomain& x, Range& y, const LocalIndexSet& localIndexSet, const Vector& coefficients, const NodeToRangeEntry& nodeToRangeEntry, ShapeFunctionValueContainer& shapeFunctionValueContainer):
        x_(x),
        y_(y),
        localIndexSet_(localIndexSet),
        coefficients_(coefficients),
        nodeToRangeEntry_(nodeToRangeEntry),
        shapeFunctionValueContainer_(shapeFunctionValueContainer)
      {}

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath treePath)
      {
        using LocalBasisRange = typename Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
        using MultiIndex = typename LocalIndexSet::MultiIndex;
        using CoefficientBlock = typename std::decay<decltype(std::declval<Vector>()[std::declval<MultiIndex>()])>::type;
        using RangeBlock = typename std::decay<decltype(nodeToRangeEntry_(node, y_))>::type;

        auto&& fe = node.finiteElement();
        auto&& localBasis = fe.localBasis();

        auto&& shapeFunctionValues = shapeFunctionValueContainer_[node];
        localBasis.evaluateFunction(x_, shapeFunctionValues);

        // Apply function space type dependent continuity preserving transformation.
        // Distinguish between H, Hdiv and Hcurl types. Use transformations:
        // H:     φ → φ
        // Hdiv:  φ → 1/|det J| J^T φ
        // Hcurl: φ → J^(-T) φ
        auto&& type = fe.functionSpaceType();
        // if type==FunctionSpace::Type::H) no further transformation needed.
        if (type==FunctionSpace::Type::Hdiv)
        {
          auto element = node.element();
          auto geometry = element.geometry();
          auto jacobianTransposed = geometry.jacobianTransposed(x_);
          auto integrationElement = geometry.integrationElement(x_);
          for (size_type i = 0; i < localBasis.size(); ++i)
          {
            auto tmp = shapeFunctionValues[i];
            jacobianTransposed.mtv(tmp, shapeFunctionValues[i]);
            shapeFunctionValues[i] /= integrationElement;
          }
        }
        else if (type==FunctionSpace::Type::Hcurl)
        {
          auto element = node.element();
          auto geometry = element.geometry();
          auto jacobianInverseTransposed = geometry.jacobianInverseTransposed(x_);
          for (size_type i = 0; i < localBasis.size(); ++i)
          {
            auto tmp = shapeFunctionValues[i];
            jacobianInverseTransposed.mtv(tmp, shapeFunctionValues[i]);
          }
        }

        // Get range entry associated to this node
        auto&& re = nodeToRangeEntry_(node, y_);

        for (size_type i = 0; i < localBasis.size(); ++i)
        {
          auto&& multiIndex = localIndexSet_.index(node.localIndex(i));

          // Get coefficient associated to i-th shape function
          auto&& c = coefficients_[multiIndex];

          // Get value of i-th shape function
          auto&& v = shapeFunctionValues[i];

          // Notice that the range entry re, the coefficient c, and the shape functions
          // value v may all be scalar, vector, matrix, or general container valued.
          // The matching of their entries is done via the multistage procedure described
          // in the class documentation of DiscreteGlobalBasisFunction.
          auto dimC = FlatVectorBackend<CoefficientBlock>::size(c);
          auto dimV = FlatVectorBackend<LocalBasisRange>::size(v);
          assert(dimC*dimV == FlatVectorBackend<RangeBlock>::size(re));
          for(size_type j=0; j<dimC; ++j)
          {
            auto&& c_j = FlatVectorBackend<CoefficientBlock>::getEntry(c, j);
            for(size_type k=0; k<dimV; ++k)
            {
              auto&& v_k = FlatVectorBackend<LocalBasisRange>::getEntry(v, k);
              FlatVectorBackend<RangeBlock>::getEntry(re, j*dimV + k) += c_j*v_k;
            }
          }
        }
      }

      const LocalDomain& x_;
      Range& y_;
      const LocalIndexSet& localIndexSet_;
      const Vector& coefficients_;
      const NodeToRangeEntry& nodeToRangeEntry_;
      ShapeFunctionValueContainer& shapeFunctionValueContainer_;
    };

    // TODO: To be tested properly!
    struct LocalJacobianEvaluateVisitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {
      LocalJacobianEvaluateVisitor(const LocalDomain& x, JacobianRange& y, const LocalIndexSet& localIndexSet, const Vector& coefficients, const NodeToRangeEntry& nodeToRangeEntry, ShapeFunctionJacobianContainer& shapeFunctionJacobianContainer):
        x_(x),
        y_(y),
        localIndexSet_(localIndexSet),
        coefficients_(coefficients),
        nodeToRangeEntry_(nodeToRangeEntry),
        shapeFunctionJacobianContainer_(shapeFunctionJacobianContainer)
      {}

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath treePath)
      {
        using LocalBasisJacobianRange = typename Node::FiniteElement::Traits::LocalBasisType::Traits::JacobianType;
        using MultiIndex = typename LocalIndexSet::MultiIndex;
        using CoefficientBlock = typename std::decay<decltype(std::declval<Vector>()[std::declval<MultiIndex>()])>::type;

        using RangeBlock = typename std::decay<decltype(nodeToRangeEntry_(node, y_))>::type;

        auto&& fe = node.finiteElement();
        auto&& localBasis = fe.localBasis();

        auto&& shapeFunctionJacobians = shapeFunctionJacobianContainer_[node];
        localBasis.evaluateJacobian(x_, shapeFunctionJacobians);

        // Apply function space type dependent continuity preserving transformation.
        // Distinguish between H, Hdiv and Hcurl types. Use transformations:
        // H:     ∇φ → J^(-T) ∇φ
        // Hdiv:  ∇φ → 1/|det J| J ∇φ J^(-1)
        // Hcurl: ∇φ → ... not implemented
        auto&& type = fe.functionSpaceType();
        auto element = node.element();
        auto geometry = element.geometry();
        if (type==FunctionSpace::Type::H)
        {
          auto jacobianInverseTransposed = geometry.jacobianInverseTransposed(x_);
          for (size_type i = 0; i < localBasis.size(); ++i)
          {
            auto referenceShapeFunctionJacobian = shapeFunctionJacobians[i];
            jacobianInverseTransposed.mv(referenceShapeFunctionJacobian[0], shapeFunctionJacobians[i][0]);
          }
        }
        else if (type==FunctionSpace::Type::Hdiv)
        {
          auto integrationElement = geometry.integrationElement(x_);
          auto jacobianInverseTransposed = geometry.jacobianInverseTransposed(x_);
          auto jacobianTransposed = geometry.jacobianTransposed(x_);
          for (size_type i = 0; i < localBasis.size(); ++i)
          {
            auto referenceShapeFunctionJacobian = shapeFunctionJacobians[i];
            jacobianTransposed.mtv(referenceShapeFunctionJacobian[0], shapeFunctionJacobians[i][0]);
            referenceShapeFunctionJacobian = shapeFunctionJacobians[i];
            jacobianInverseTransposed.mtv(referenceShapeFunctionJacobian[0], shapeFunctionJacobians[i][0]);
            shapeFunctionJacobians[i] /= integrationElement;
          }
        }
        else if (type==FunctionSpace::Type::Hcurl)
        {
          DUNE_THROW(Dune::NotImplemented, "Derivatives for curl function not supported by dune-functions.");
        }

        // Get range entry associated to this node
        auto&& re = nodeToRangeEntry_(node, y_);

        for (size_type i = 0; i < localBasis.size(); ++i)
        {
          auto&& multiIndex = localIndexSet_.index(node.localIndex(i));

          // Get coefficient associated to i-th shape function
          auto&& c = coefficients_[multiIndex];

          // Get value of i-th shape function
          auto&& v = shapeFunctionJacobians[i];

          // Notice that the range entry re, the coefficient c, and the shape functions
          // value v may all be scalar, vector, matrix, or general container valued.
          // The matching of their entries is done via the multistage procedure described
          // in the class documentation of DiscreteGlobalBasisFunction.
          auto dimC = FlatVectorBackend<CoefficientBlock>::size(c);
          auto dimV = FlatVectorBackend<LocalBasisJacobianRange>::size(v);
          assert(dimC*dimV == FlatVectorBackend<RangeBlock>::size(re));
          for(size_type j=0; j<dimC; ++j)
          {
            auto&& c_j = FlatVectorBackend<CoefficientBlock>::getEntry(c, j);
            for(size_type k=0; k<dimV; ++k)
            {
              auto&& v_k = FlatVectorBackend<LocalBasisJacobianRange>::getEntry(v, k);
              FlatVectorBackend<RangeBlock>::getEntry(re, j*dimV + k) += c_j*v_k;
            }
          }
        }
      }

      const LocalDomain& x_;
      JacobianRange& y_;
      const LocalIndexSet& localIndexSet_;
      const Vector& coefficients_;
      const NodeToRangeEntry& nodeToRangeEntry_;
      ShapeFunctionJacobianContainer& shapeFunctionJacobianContainer_;
    };

    struct LocalDofVisitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {
      LocalDofVisitor (std::vector<RangeField>& localDofs, const LocalIndexSet& localIndexSet, const Vector& coefficients) :
        localDofs_(localDofs),
        localIndexSet_(localIndexSet),
        coefficients_(coefficients),
        counter(0)
      {}

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath treePath)
      {
        using MultiIndex = typename LocalIndexSet::MultiIndex;
        using CoefficientBlock = typename std::decay<decltype(std::declval<Vector>()[std::declval<MultiIndex>()])>::type;

        auto&& fe = node.finiteElement();
        auto&& localBasis = fe.localBasis();

        for (size_type i = 0; i < localBasis.size(); ++i)
        {
          auto&& multiIndex = localIndexSet_.index(node.localIndex(i));

          // Get coefficient associated to i-th shape function
          auto&& c = coefficients_[multiIndex];

          // Write modus depends on index type (interleafed, ...)
          localDofs_[counter++] = c;

          auto dimC = FlatVectorBackend<CoefficientBlock>::size(c);
          if (dimC != 1 )std::cout << "[Warning] getLocalDofs() does only work for scalar local coefficients. (dune-functions)" << std::endl;
        }
      }

      std::vector<RangeField>& localDofs_;
      const LocalIndexSet& localIndexSet_;
      const Vector& coefficients_;
      size_t counter=0;
    };

  public:

    using GlobalFunction = DiscreteGlobalBasisFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using JacobianRange = GlobalFunction::JacobianRange;
    using Element = GlobalFunction::Element;

    LocalFunction(const DiscreteGlobalBasisFunction& globalFunction)
      : globalFunction_(&globalFunction)
      , localBasisView_(globalFunction.basis().localView())
      , localIndexSet_(globalFunction.basis().localIndexSet())
    {
      // Here we assume that the tree can be accessed, traversed,
      // and queried for tree indices even in unbound state.
      subTree_ = &TypeTree::child(localBasisView_.tree(), globalFunction_->treePath());
      shapeFunctionValueContainer_.init(*subTree_);
      shapeFunctionJacobianContainer_.init(*subTree_);
//      localDoFs_.reserve(localBasisView_.maxSize());
    }

    LocalFunction(const LocalFunction& other)
      : globalFunction_(other.globalFunction_)
      , localBasisView_(globalFunction_->basis().localView())
      , localIndexSet_(globalFunction_->basis().localIndexSet())
    {
      // Here we assume that the tree can be accessed, traversed,
      // and queried for tree indices even in unbound state.
      subTree_ = &TypeTree::child(localBasisView_.tree(), globalFunction_->treePath());
      shapeFunctionValueContainer_.init(*subTree_);
      shapeFunctionJacobianContainer_.init(*subTree_);
    }

    LocalFunction operator=(const LocalFunction& other)
    {
      globalFunction_ = other.globalFunction_;
      localBasisView_ = other.localBasisView_;
      localIndexSet_ = other.localIndexSet_;
      subTree_ = &TypeTree::child(localBasisView_.tree(), globalFunction_->treePath());

      // Here we assume that the tree can be accessed, traversed,
      // and queried for tree indices even in unbound state.
      shapeFunctionValueContainer_.init(*subTree_);
      shapeFunctionJacobianContainer_.init(*subTree_);
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluate()
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element)
    {
      localBasisView_.bind(element);
      localIndexSet_.bind(localBasisView_);

      // Read dofs associated to bound element
//      localDoFs_.resize(subTree_->size());
//      for (size_type i = 0; i < subTree_->size(); ++i)
//        localDoFs_[i] = globalFunction_->dofs()[localIndexSet_.index(i)];
    }

    void unbind()
    {
      localIndexSet_.unbind();
      localBasisView_.unbind();
    }

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

      LocalEvaluateVisitor localEvaluateVisitor(x, y, localIndexSet_, globalFunction_->dofs(), globalFunction_->nodeToRangeEntry(), shapeFunctionValueContainer_);
      TypeTree::applyToTree(*subTree_, localEvaluateVisitor);

      return y;
    }

    /**
     * \brief Evaluate derivative of LocalFunction at bound element.
     */
    JacobianRange derivative(const Domain& x) const
    {
      JacobianRange y;
      for (size_t j=0; j<y.size(); ++j)
        y[j] = 0.0;

      LocalJacobianEvaluateVisitor localJacobianEvaluateVisitor(x, y, localIndexSet_, globalFunction_->dofs(), globalFunction_->nodeToRangeEntry(), shapeFunctionJacobianContainer_);
      TypeTree::applyToTree(*subTree_, localJacobianEvaluateVisitor);

      return y;
    }

    std::vector<RangeField> localDofs() const
    {
      std::vector<RangeField> localDofs(subTree_->size());

      LocalDofVisitor localDofVisitor(localDofs, localIndexSet_, globalFunction_->dofs());
      TypeTree::applyToTree(*subTree_, localDofVisitor);

      return localDofs;
    }

    const Element& localContext() const
    {
      return localBasisView_.element();
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:

    const DiscreteGlobalBasisFunction* globalFunction_;
    LocalBasisView localBasisView_;
    LocalIndexSet localIndexSet_;

    mutable ShapeFunctionValueContainer shapeFunctionValueContainer_;
    mutable ShapeFunctionJacobianContainer shapeFunctionJacobianContainer_;
//    std::vector<typename V::value_type> localDoFs_;
    const SubTree* subTree_;
  };

  DiscreteGlobalBasisFunction(const Basis & basis, const TreePath& treePath, const V & coefficients, const NodeToRangeEntry& nodeToRangeEntry) :
    entitySet_(basis.gridView()),
    basis_(stackobject_to_shared_ptr(basis)),
    treePath_(treePath),
    coefficients_(stackobject_to_shared_ptr(coefficients)),
    nodeToRangeEntry_(stackobject_to_shared_ptr(nodeToRangeEntry))
  {}

  DiscreteGlobalBasisFunction(std::shared_ptr<const Basis> basis, const TreePath& treePath, std::shared_ptr<const V> coefficients, std::shared_ptr<const NodeToRangeEntry> nodeToRangeEntry) :
    entitySet_(basis->gridView()),
    basis_(basis),
    treePath_(treePath),
    coefficients_(coefficients),
    nodeToRangeEntry_(nodeToRangeEntry)
  {}

  const Basis& basis() const
  {
    return *basis_;
  }

  const TreePath& treePath() const
  {
    return treePath_;
  }

  const V& dofs() const
  {
    return *coefficients_;
  }

  const NodeToRangeEntry& nodeToRangeEntry() const
  {
    return *nodeToRangeEntry_;
  }

  // TODO: Implement this using hierarchic search
  Range operator() (const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend typename Traits::DerivativeInterface derivative(const DiscreteGlobalBasisFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend LocalFunction localFunction(const DiscreteGlobalBasisFunction& t)
  {
    return LocalFunction(t);
  }

  /**
   * \brief Get associated EntitySet
   */
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:

  EntitySet entitySet_;
  std::shared_ptr<const Basis> basis_;
  const TreePath treePath_;
  std::shared_ptr<const V> coefficients_;
  std::shared_ptr<const NodeToRangeEntry> nodeToRangeEntry_;
};



template<typename R, typename B, typename TP, typename V>
auto makeDiscreteGlobalBasisFunction(B&& basis, const TP& treePath, V&& vector)
{
  using Basis = std::decay_t<B>;
  using Vector = std::decay_t<V>;
  using NTREM = DefaultNodeToRangeMap<typename TypeTree::ChildForTreePath<typename Basis::LocalView::Tree, TP>>;
  auto nodeToRangeEntryPtr = std::make_shared<NTREM>(makeDefaultNodeToRangeMap(basis, treePath));
  auto basisPtr = Dune::wrap_or_move(std::forward<B>(basis));
  auto vectorPtr = Dune::wrap_or_move(std::forward<V>(vector));
  return DiscreteGlobalBasisFunction<Basis, TP, Vector, NTREM, R>(basisPtr, treePath, vectorPtr, nodeToRangeEntryPtr);
}



template<typename R, typename B, typename V>
auto makeDiscreteGlobalBasisFunction(B&& basis, V&& vector)
{
  return makeDiscreteGlobalBasisFunction<R>(std::forward<B>(basis), TypeTree::hybridTreePath(), std::forward<V>(vector));
}



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
