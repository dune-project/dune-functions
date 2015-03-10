// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH

#include <memory>

#include <dune/common/shared_ptr.hh>

#include <dune/functions/gridfunctions/gridviewentityset.hh>
//#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {
namespace Functions {

template<typename Basis, typename V>
class DiscreteScalarGlobalBasisFunction
{
public:
  using GridView = typename Basis::GridView;
  using EntitySet = GridViewEntitySet<GridView, 0>;

  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeType;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  class LocalFunction
  {
    using LocalBasisView = typename Basis::LocalView;
    using LocalIndexSet = typename Basis::IndexSet::LocalIndexSet;
    using size_type = typename LocalBasisView::Tree::size_type;

  public:

    using GlobalFunction = DiscreteScalarGlobalBasisFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using Element = GlobalFunction::Element;

    LocalFunction(const DiscreteScalarGlobalBasisFunction& globalFunction)
      : globalFunction_(globalFunction)
      , localBasisView_(globalFunction.basis().localView())
      , localIndexSet_(globalFunction.indexSet_.localIndexSet())
    {
      localDoFs_.reserve(localBasisView_.maxSize());
    }

//    typename EBase::DerivativeBasePointer derivative() const
//    {
//      DUNE_THROW(NotImplemented,"derivative not implemented");
//    }

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

      auto& tree = localBasisView_.tree();

      // Read dofs associated to bound element
      localDoFs_.resize(localIndexSet_.size());
      for (size_type i = 0; i < localIndexSet_.size(); ++i)
        localDoFs_[i] = globalFunction_.dofs()[localIndexSet_.index(i)[0]];
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
     * you have to call bind() again in order to make evaluate()
     * usable.
     */
    void evaluate(const Domain& coord, Range& r) const
    {
      std::vector<Range> shapeFunctionValues;
      auto& basis = localBasisView_.tree().finiteElement().localBasis();
      basis.evaluateFunction(coord,shapeFunctionValues);
      r = 0;
      for (size_type i = 0; i < basis.size(); ++i)
        r += localDoFs_[i] * shapeFunctionValues[i];
    }

    Range operator()(const Domain& x) const
    {
      Range y;
      evaluate(x,y);
      return y;
    }

    const Element& localContext() const
    {
      return localBasisView_.element();
    }

  private:

    const DiscreteScalarGlobalBasisFunction& globalFunction_;
    LocalBasisView localBasisView_;
    LocalIndexSet localIndexSet_;
    std::vector<typename V::value_type> localDoFs_;
  };

  DiscreteScalarGlobalBasisFunction(const Basis & basis, const V & dofs)
    : basis_(stackobject_to_shared_ptr(basis))
    , dofs_(stackobject_to_shared_ptr(dofs))
    , indexSet_(basis.indexSet())
  {}

  DiscreteScalarGlobalBasisFunction(std::shared_ptr<Basis> basis, std::shared_ptr<V> dofs)
    : basis_(basis)
    , dofs_(dofs)
    , indexSet_(basis.indexSet())
  {}

  const Basis& basis() const
  {
    return *basis_;
  }

  const V& dofs() const
  {
    return *dofs_;
  }

  // TODO: Implement this using hierarchic search
  Range operator() (const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend LocalFunction localFunction(const DiscreteScalarGlobalBasisFunction& t)
  {
    return LocalFunction(t);
  }

private:

  std::shared_ptr<const Basis> basis_;
  std::shared_ptr<const V> dofs_;
  typename Basis::IndexSet indexSet_;

};

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH
