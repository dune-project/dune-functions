// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH

#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {
namespace Functions {

template<typename Basis, typename V>
class DiscreteScalarGlobalBasisFunction
  : public GridViewFunction<typename Basis::GridView,
                            typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeType
                            >
{

public:

  typedef GridViewFunction<
    typename Basis::GridView,
    typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeType
    > Base;

  typedef typename Base::Element Element;

  class LocalFunction
    : public Base::LocalFunction
  {

    typedef typename Base::LocalFunction EBase;
    typedef typename Basis::LocalView LocalBasisView;
    typedef typename LocalBasisView::Tree::size_type size_type;

  public:

    typedef typename EBase::LocalContext Element;
    typedef typename EBase::Domain Domain;
    typedef typename EBase::Range Range;

    LocalFunction(const DiscreteScalarGlobalBasisFunction& globalFunction)
      : globalFunction_(globalFunction)
      , localBasisView_(globalFunction.basis().localView())
    {
      indexCache_.reserve(globalFunction_.basis().maxLocalSize());
      dofCache_.reserve(globalFunction_.basis().maxLocalSize());
    }

    virtual typename EBase::DerivativeBasePointer derivative() const DUNE_FINAL
    {
      DUNE_THROW(NotImplemented,"derivative not implemented");
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluate()
     * and after changes to the coefficient vector.
     */
    virtual void bind(const Element& element) DUNE_FINAL
    {
      localBasisView_.bind(element);

      auto& tree = localBasisView_.tree();

      // Stricly speaking we do not need to cache the indices if we
      // already cache the dofs.
      indexCache_.resize(tree.localSize());
      tree.generateMultiIndices(indexCache_.begin());

      // Cache the dofs associated to bound element
      dofCache_.resize(tree.localSize());
      for (size_type i = 0; i < tree.localSize(); ++i)
        dofCache_[i] = globalFunction_.dofs()[indexCache_[i][0]];
    }

    virtual void unbind() DUNE_FINAL
    {
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
    virtual void evaluate(const Domain& coord, Range& r) const DUNE_FINAL
    {
      std::vector<Range> shapeFunctionValues;
      auto& basis = localBasisView_.tree().finiteElement().localBasis();
      basis.evaluateFunction(coord,shapeFunctionValues);
      r = 0;
      for (size_type i = 0; i < basis.size(); ++i)
        r += dofCache_[i] * shapeFunctionValues[i];
    }

    virtual const Element& localContext() const DUNE_FINAL
    {
      return localBasisView_.element();
    }

  private:

    const DiscreteScalarGlobalBasisFunction& globalFunction_;
    LocalBasisView localBasisView_;
    std::vector<typename LocalBasisView::MultiIndex> indexCache_;
    std::vector<typename V::value_type> dofCache_;

  };

  DiscreteScalarGlobalBasisFunction(std::shared_ptr<Basis> basis, std::shared_ptr<V> dofs)
    : Base(basis->gridView())
    , basis_(basis)
    , dofs_(dofs)
  {}

  virtual typename Base::LocalFunctionBasePointer localFunction() const DUNE_FINAL
  {
    return std::make_shared<LocalFunction>(*this);
  }

  const Basis& basis() const
  {
    return *basis_;
  }

  const V& dofs() const
  {
    return *dofs_;
  }

  virtual typename Base::DerivativeBasePointer derivative() const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"derivative not implemented yet");
  }

  // TODO: Implement this using hierarchic search
  virtual void evaluate(const typename Base::Domain& domain, typename Base::Range& r) const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }



private:

  std::shared_ptr<Basis> basis_;
  std::shared_ptr<V> dofs_;

};

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH
