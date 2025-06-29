// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH

#include <memory>
#include <optional>

#include <dune/common/typetraits.hh>

#include <dune/grid/utility/hierarchicsearch.hh>

#include <dune/typetree/treecontainer.hh>

#include <dune/functions/functionspacebases/hierarchicnodetorangemap.hh>
#include <dune/functions/functionspacebases/flatvectorview.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/backends/concepts.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

namespace Dune {
namespace Functions {


namespace ImplDoc {

template<typename B, typename V, typename NTRE>
class DiscreteGlobalBasisFunctionBase
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

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

protected:

  // This collects all data that is shared by all related
  // global and local functions. This way we don't need to
  // keep track of it individually.
  struct Data
  {
    EntitySet entitySet;
    std::shared_ptr<const Basis> basis;
    std::shared_ptr<const Vector> coefficients;
    std::shared_ptr<const NodeToRangeEntry> nodeToRangeEntry;
  };

public:
  class LocalFunctionBase
  {
    using LocalView = typename Basis::LocalView;
    using size_type = typename Tree::size_type;

  public:
    using Domain = LocalDomain;
    using Element = typename EntitySet::Element;

  protected:
    LocalFunctionBase(const std::shared_ptr<const Data>& data)
      : data_(data)
      , localView_(data_->basis->localView())
    {
      localDoFs_.reserve(localView_.maxSize());
    }

    /**
     * \brief Copy-construct the local-function.
     *
     * This copy-constructor copies the cached local DOFs only
     * if the `other` local-function is bound to an element.
     **/
    LocalFunctionBase(const LocalFunctionBase& other)
      : data_(other.data_)
      , localView_(other.localView_)
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
    LocalFunctionBase& operator=(const LocalFunctionBase& other)
    {
      data_ = other.data_;
      localView_ = other.localView_;
      if (bound())
        localDoFs_ = other.localDoFs_;
      return *this;
    }

  public:
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
      const auto& dofs = *data_->coefficients;
      for (size_type i = 0; i < localView_.tree().size(); ++i)
      {
        // For a subspace basis the index-within-tree i
        // is not the same as the localIndex within the
        // full local view.
        size_t localIndex = localView_.tree().localIndex(i);
        localDoFs_[localIndex] = dofs[localView_.index(localIndex)];
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

    //! Return the element the local-function is bound to.
    const Element& localContext() const
    {
      return localView_.element();
    }

  protected:

    template<class To, class From>
    void assignWith(To& to, const From& from) const
    {
        auto from_flat = flatVectorView(from);
        auto to_flat = flatVectorView(to);
        assert(from_flat.size() == to_flat.size());
        for (size_type i = 0; i < to_flat.size(); ++i)
          to_flat[i] = from_flat[i];
    }

    template<class Node, class TreePath, class Range>
    decltype(auto) nodeToRangeEntry(const Node& node, const TreePath& treePath, Range& y) const
    {
      return (*data_->nodeToRangeEntry)(node, treePath, y);
    }

    std::shared_ptr<const Data> data_;
    LocalView localView_;
    std::vector<Coefficient> localDoFs_;
  };

protected:
  DiscreteGlobalBasisFunctionBase(const std::shared_ptr<const Data>& data)
    : data_(data)
  {
    /* Nothing. */
  }

public:

  //! Return a const reference to the stored basis.
  const Basis& basis() const
  {
    return *data_->basis;
  }

  //! Return the coefficients of this discrete function by reference.
  const Vector& dofs() const
  {
    return *data_->coefficients;
  }

  //! Return the stored node-to-range map.
  const NodeToRangeEntry& nodeToRangeEntry() const
  {
    return *data_->nodeToRangeEntry;
  }

  //! Get associated set of entities the local-function can be bound to.
  const EntitySet& entitySet() const
  {
    return data_->entitySet;
  }

protected:
  std::shared_ptr<const Data> data_;
};

} // namespace ImplDoc



template<typename DGBF>
class DiscreteGlobalBasisFunctionDerivative;

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
 * 1.Each leaf node in the local ansatz subtree is associated to an
 *   entry `RE` of the range-type via the given node-to-range-entry-map.
 *   Based on this mapping each node is processed independently in the
 *   following way:
 *
 * 2.Now let the coefficients type `C` per basis function be `dim_C`-dimensional.
 *   Then we compute the dim(C) linear combinations (one for each coefficient
 *   index) of the shape function values with type `V` independently storing them in
 *   a `std::array<V,dim_C>`.
 *
 * 3.Finally the resulting array of function values is assigned to the
 *   nodal range entry `RE`. Since both types may be different their entries
 *   are mapped to one another via `flatVectorView()`. This will recursive
 *   enumerate the entries of the types in lexicographic order (unless
 *   `flatVectorView` is specialized differently for a certain type).
 *
 * As a consequence the nodal range entry is required to have a total
 * dimension `dim_RE = dim_C * dim_V` and to be compatible with `flatVectorView()`.
 *
 * \tparam B Type of global basis
 * \tparam V Type of coefficient vectors
 * \tparam NTRE Type of node-to-range-entry-map that associates each leaf node in the local ansatz subtree with an entry in the range type
 * \tparam R Range type of this function
 *
 * \b Requirements:
 * - Type `R` is default constructible
 * - default constructed `R` has correct "shape" to hold the components of the range and
 *   to be accessible using the node-to-range-entry map `NTRE`.
 * - This requirements applies recursively to the derivative ranges.
 */
template<typename B, typename V,
  typename NTRE = HierarchicNodeToRangeMap,
  typename R = typename V::value_type>
class DiscreteGlobalBasisFunction
  : public ImplDoc::DiscreteGlobalBasisFunctionBase<B, V, NTRE>
{
  using Base = ImplDoc::DiscreteGlobalBasisFunctionBase<B, V, NTRE>;
  using Data = typename Base::Data;

public:
  using Basis = typename Base::Basis;
  using Vector = typename Base::Vector;

  using Domain = typename Base::Domain;
  using Range = R;

  using Traits = Imp::GridFunctionTraits<Range(Domain), typename Base::EntitySet, DefaultDerivativeTraits, 16>;

private:

  template<class Node>
  using LocalBasisRange = typename Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
  template<class Node>
  using NodeData = typename std::vector<LocalBasisRange<Node>>;
  using PerNodeEvaluationBuffer = typename TypeTree::TreeContainer<NodeData, typename Base::Tree>;

public:
  class LocalFunction
    : public Base::LocalFunctionBase
  {
    using LocalBase = typename Base::LocalFunctionBase;
    using size_type = typename Base::Tree::size_type;
    using LocalBase::nodeToRangeEntry;

  public:

    using GlobalFunction = DiscreteGlobalBasisFunction;
    using Domain = typename LocalBase::Domain;
    using Range = GlobalFunction::Range;
    using Element = typename LocalBase::Element;

    //! Create a local-function from the associated grid-function
    LocalFunction(const DiscreteGlobalBasisFunction& globalFunction)
      : LocalBase(globalFunction.data_)
      , evaluationBuffer_(this->localView_.tree())
    {
      /* Nothing. */
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

      TypeTree::forEachLeafNode(this->localView_.tree(), [&](auto&& node, auto&& treePath) {
        if (node.empty())
          return;
        const auto& fe = node.finiteElement();
        const auto& localBasis = fe.localBasis();
        auto& shapeFunctionValues = evaluationBuffer_[treePath];

        localBasis.evaluateFunction(x, shapeFunctionValues);

        // Compute linear combinations of basis function jacobian.
        // Non-scalar coefficients of dimension coeffDim are handled by
        // processing the coeffDim linear combinations independently
        // and storing them as entries of an array.
        using Value = LocalBasisRange< std::decay_t<decltype(node)> >;
        static constexpr auto coeffDim = decltype(flatVectorView(this->localDoFs_[node.localIndex(0)]).size())::value;
        auto values = std::array<Value, coeffDim>{};
        istlVectorBackend(values) = 0;
        for (size_type i = 0; i < localBasis.size(); ++i)
        {
          auto c = flatVectorView(this->localDoFs_[node.localIndex(i)]);
          for (std::size_t j = 0; j < coeffDim; ++j)
            values[j].axpy(c[j], shapeFunctionValues[i]);
        }

        // Assign computed values to node entry of range.
        // Types are matched using the lexicographic ordering provided by flatVectorView.
        LocalBase::assignWith(nodeToRangeEntry(node, treePath, y), values);
      });

      return y;
    }

    //! Local function of the derivative
    friend typename DiscreteGlobalBasisFunctionDerivative<DiscreteGlobalBasisFunction>::LocalFunction derivative(const LocalFunction& lf)
    {
      auto dlf = localFunction(DiscreteGlobalBasisFunctionDerivative<DiscreteGlobalBasisFunction>(lf.data_));
      if (lf.bound())
        dlf.bind(lf.localContext());
      return dlf;
    }

  private:
    mutable PerNodeEvaluationBuffer evaluationBuffer_;
  };

  //! Create a grid-function, by wrapping the arguments in `std::shared_ptr`.
  template<class B_T, class V_T, class NTRE_T>
  DiscreteGlobalBasisFunction(B_T && basis, V_T && coefficients, NTRE_T&& nodeToRangeEntry)
    : Base(std::make_shared<Data>(Data{{basis.gridView()}, wrap_or_move(std::forward<B_T>(basis)), wrap_or_move(std::forward<V_T>(coefficients)), wrap_or_move(std::forward<NTRE_T>(nodeToRangeEntry))}))
  {}

  //! Create a grid-function, by moving the arguments in `std::shared_ptr`.
  DiscreteGlobalBasisFunction(std::shared_ptr<const Basis> basis, std::shared_ptr<const V> coefficients, std::shared_ptr<const typename Base::NodeToRangeEntry> nodeToRangeEntry)
    : Base(std::make_shared<Data>(Data{{basis->gridView()}, basis, coefficients, nodeToRangeEntry}))
  {}

  /** \brief Evaluate at a point given in world coordinates
   *
   * \warning This has to find the element that the evaluation point is in.
   *   It is therefore very slow.
   */
  Range operator() (const Domain& x) const
  {
    HierarchicSearch search(this->data_->basis->gridView().grid(), this->data_->basis->gridView().indexSet());

    const auto e = search.findEntity(x);
    auto localThis = localFunction(*this);
    localThis.bind(e);
    return localThis(e.geometry().local(x));
  }

  //! Derivative of the `DiscreteGlobalBasisFunction`
  friend DiscreteGlobalBasisFunctionDerivative<DiscreteGlobalBasisFunction> derivative(const DiscreteGlobalBasisFunction& f)
  {
    return DiscreteGlobalBasisFunctionDerivative<DiscreteGlobalBasisFunction>(f.data_);
  }

  /**
   * \brief Construct local function from a DiscreteGlobalBasisFunction.
   *
   * The obtained local function satisfies the concept
   * `Dune::Functions::Concept::LocalFunction`. It must be bound
   * to an entity from the entity set of the DiscreteGlobalBasisFunction
   * before it can be used.
   */
  friend LocalFunction localFunction(const DiscreteGlobalBasisFunction& t)
  {
    return LocalFunction(t);
  }
};


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
 *            for more details. The type `R` must be default constructible and
 *            by this default construction provides a well resized container to be
 *            filled in the evaluation method of the grid function.
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


/**
 * \brief Derivative of a `DiscreteGlobalBasisFunction`
 *
 * Function returning the derivative of the given `DiscreteGlobalBasisFunction`
 * with respect to global coordinates.
 *
 * The function handles the mapping of coefficient blocks and basis function values
 * to range entries analogously to the `DiscreteGlobalBasisFunction`. This mapping
 * is implemented with the same algorithm but with values replace by jacobian.
 *
 * \ingroup FunctionImplementations
 *
 * \tparam DGBF instance of the `DiscreteGlobalBasisFunction` this is a derivative of
 */
template<typename DGBF>
class DiscreteGlobalBasisFunctionDerivative
  : public ImplDoc::DiscreteGlobalBasisFunctionBase<typename DGBF::Basis, typename DGBF::Vector, typename DGBF::NodeToRangeEntry>
{
  using Base = ImplDoc::DiscreteGlobalBasisFunctionBase<typename DGBF::Basis, typename DGBF::Vector, typename DGBF::NodeToRangeEntry>;
  using Data = typename Base::Data;

public:
  using DiscreteGlobalBasisFunction = DGBF;

  using Basis = typename Base::Basis;
  using Vector = typename Base::Vector;

  using Domain = typename Base::Domain;
  using Range = typename SignatureTraits<typename DiscreteGlobalBasisFunction::Traits::DerivativeInterface>::Range;

  using Traits = Imp::GridFunctionTraits<Range(Domain), typename Base::EntitySet, DefaultDerivativeTraits, 16>;

private:

  template<class Node>
  using LocalBasisRange = typename Node::FiniteElement::Traits::LocalBasisType::Traits::JacobianType;
  template<class Node>
  using NodeData = typename std::vector< LocalBasisRange<Node> >;
  using PerNodeEvaluationBuffer = typename TypeTree::TreeContainer<NodeData, typename Base::Tree>;

public:

  /**
   * \brief local function evaluating the derivative in reference coordinates
   *
   * Note that the function returns the derivative with respect to global
   * coordinates even when the point is given in reference coordinates on
   * an element.
   */
  class LocalFunction
    : public Base::LocalFunctionBase
  {
    using LocalBase = typename Base::LocalFunctionBase;
    using size_type = typename Base::Tree::size_type;
    using LocalBase::nodeToRangeEntry;

  public:
    using GlobalFunction = DiscreteGlobalBasisFunctionDerivative;
    using Domain = typename LocalBase::Domain;
    using Range = GlobalFunction::Range;
    using Element = typename LocalBase::Element;

    //! Create a local function from the associated grid function
    LocalFunction(const GlobalFunction& globalFunction)
      : LocalBase(globalFunction.data_)
      , evaluationBuffer_(this->localView_.tree())
    {
      /* Nothing. */
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before `operator()`
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element)
    {
      LocalBase::bind(element);
      geometry_.emplace(element.geometry());
    }

    //! Unbind the local-function.
    void unbind()
    {
      geometry_.reset();
      LocalBase::unbind();
    }

    /**
     * \brief Evaluate this local-function in coordinates `x` in the bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     *
     * Note that the function returns the derivative with respect to global
     * coordinates even when the point is given in reference coordinates on
     * an element.
     */
    Range operator()(const Domain& x) const
    {
      Range y;
      istlVectorBackend(y) = 0;

      const auto& jacobianInverse = geometry_->jacobianInverse(x);

      TypeTree::forEachLeafNode(this->localView_.tree(), [&](auto&& node, auto&& treePath) {
        if (node.empty())
          return;
        const auto& fe = node.finiteElement();
        const auto& localBasis = fe.localBasis();
        auto& shapeFunctionJacobians = evaluationBuffer_[treePath];

        localBasis.evaluateJacobian(x, shapeFunctionJacobians);

        // Compute linear combinations of basis function jacobian.
        // Non-scalar coefficients of dimension coeffDim are handled by
        // processing the coeffDim linear combinations independently
        // and storing them as entries of an array.
        using RefJacobian = LocalBasisRange< std::decay_t<decltype(node)> >;
        static constexpr auto coeffDim = decltype(flatVectorView(this->localDoFs_[node.localIndex(0)]).size())::value;
        auto refJacobians = std::array<RefJacobian, coeffDim>{};
        istlVectorBackend(refJacobians) = 0;
        for (size_type i = 0; i < localBasis.size(); ++i)
        {
          auto c = flatVectorView(this->localDoFs_[node.localIndex(i)]);
          for (std::size_t j = 0; j < coeffDim; ++j)
            refJacobians[j].axpy(c[j], shapeFunctionJacobians[i]);
        }

        // Transform Jacobians form local to global coordinates.
        using Jacobian = decltype(refJacobians[0] * jacobianInverse);
        auto jacobians = std::array<Jacobian, coeffDim>{};
        std::transform(
          refJacobians.begin(), refJacobians.end(), jacobians.begin(),
          [&](const auto& refJacobian) { return refJacobian * jacobianInverse; });

        // Assign computed Jacobians to node entry of range.
        // Types are matched using the lexicographic ordering provided by flatVectorView.
        LocalBase::assignWith(nodeToRangeEntry(node, treePath, y), jacobians);
      });

      return y;
    }

    //! Not implemented
    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction&)
    {
      DUNE_THROW(NotImplemented, "derivative of derivative is not implemented");
    }

  private:
    mutable PerNodeEvaluationBuffer evaluationBuffer_;
    std::optional<typename Element::Geometry> geometry_;
  };

  /**
   * \brief create object from `DiscreateGlobalBasisFunction` data
   *
   * Please call `derivative(discreteGlobalBasisFunction)` to create an instance
   * of this class.
   */
  DiscreteGlobalBasisFunctionDerivative(const std::shared_ptr<const Data>& data)
    : Base(data)
  {
    /* Nothing. */
  }

  /** \brief Evaluate the discrete grid-function derivative in global coordinates
   *
   * \warning This has to find the element that the evaluation point is in.
   *   It is therefore very slow.
   */
  Range operator()(const Domain& x) const
  {
    HierarchicSearch search(this->data_->basis->gridView().grid(), this->data_->basis->gridView().indexSet());

    const auto e = search.findEntity(x);
    auto localThis = localFunction(*this);
    localThis.bind(e);
    return localThis(e.geometry().local(x));
  }

  friend typename Traits::DerivativeInterface derivative(const DiscreteGlobalBasisFunctionDerivative& f)
  {
    DUNE_THROW(NotImplemented, "derivative of derivative is not implemented");
  }

  //! Construct local function from a `DiscreteGlobalBasisFunctionDerivative`
  friend LocalFunction localFunction(const DiscreteGlobalBasisFunctionDerivative& f)
  {
    return LocalFunction(f);
  }
};


} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETEGLOBALBASISFUNCTIONS_HH
