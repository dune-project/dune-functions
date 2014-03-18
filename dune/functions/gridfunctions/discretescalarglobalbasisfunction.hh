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
      : _globalFunction(globalFunction)
      , _localBasisView(globalFunction.basis().localView())
    {}

    virtual typename EBase::DerivativeBasePointer derivative() const DUNE_FINAL
    {
      DUNE_THROW(NotImplemented,"derivative not implemented");
    }

    virtual void bind(const Element& element) DUNE_FINAL
    {
      _localBasisView.bind(element);
    }

    virtual void unbind() DUNE_FINAL
    {
      _localBasisView.unbind();
    }

    virtual void evaluate(const Domain& coord, Range& r) const DUNE_FINAL
    {
      std::vector<Range> shapeFunctionValues;
      auto& basis = _localBasisView.tree().finiteElement().localBasis();
      basis.evaluateFunction(coord,shapeFunctionValues);
      std::vector<typename LocalBasisView::MultiIndex> globalIndices(basis.size());
      _localBasisView.tree().generateMultiIndices(globalIndices.begin());
      r = 0;
      for (size_type i = 0; i < basis.size(); ++i)
        r += _globalFunction.dofs()[globalIndices[i][0]] * shapeFunctionValues[i];
    }

    virtual const Element& localContext() const DUNE_FINAL
    {
      return _localBasisView.element();
    }

  private:

    const DiscreteScalarGlobalBasisFunction& _globalFunction;
    LocalBasisView _localBasisView;

  };

  DiscreteScalarGlobalBasisFunction(std::shared_ptr<Basis> basis, std::shared_ptr<V> dofs)
    : Base(basis->gridView())
    , _basis(basis)
    , _dofs(dofs)
  {}

  virtual typename Base::LocalFunctionBasePointer localFunction() const DUNE_FINAL
  {
    return std::make_shared<LocalFunction>(*this);
  }

  const Basis& basis() const
  {
    return *_basis;
  }

  const V& dofs() const
  {
    return *_dofs;
  }

  virtual typename Base::DerivativeBasePointer derivative() const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"derivative not implemented yet");
  }

  virtual void evaluate(const typename Base::Domain& domain, typename Base::Range& r) const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }



private:

  std::shared_ptr<Basis> _basis;
  std::shared_ptr<V> _dofs;

};

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH
