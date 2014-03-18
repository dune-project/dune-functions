// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1FUNCTIONSPACEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1FUNCTIONSPACEBASIS_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>


namespace Dune {
namespace Functions {


template<typename GV>
class PQ1FunctionSpaceBasisLocalView;

template<typename GV>
class PQ1FunctionSpaceBasisLeafNode;

//template<typename T>
//struct FunctionSpaceBasisTraits;
//{
//  typedef calc<T>::LocalBasisView;
//};


template<typename GV>
class PQ1FunctionSpaceBasis
{

public:
  typedef GV GridView;
  typedef std::size_t size_type;
  typedef PQ1FunctionSpaceBasisLocalView<GV> LocalView;

  static const int dim = GV::Grid::dimension;

//  typedef Dune::ReservedVector<size_type, 1> MultiIndex;
  typedef std::array<size_type, 1> MultiIndex;

  PQ1FunctionSpaceBasis(const GridView& gv) :
    gridView_(gv)
  {}

  /**
   * \brief
   */
  const GridView& gridView() const
  {
    return gridView_;
  }


  /**
   * \brief maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   * max{GridViewLocalBasisView(e).tree().size() | e in GridView}
   */
  size_type maxLocalSize() const
  {
    return 1<<dim;
  }

  //! Return number possible values for next position in empty multi index
  size_type subIndexCount() const
  {
    return gridView_.size(dim);
  }

  //! Return number possible values for next position in multi index
  size_type subIndexCount(const MultiIndex& index) const
  {
    return gridView_.size(dim);
  }

  /**
   * \brief Return local view for basis
   *
   * Perhaps we must move the construction outside
   * of the global basis in order to calm the compiler
   * when instanciating the TMP constructing the local view.
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

protected:
  const GridView gridView_;
};


template<typename GV>
class PQ1FunctionSpaceBasisLocalView
{
public:

  typedef PQ1FunctionSpaceBasis<GV> GlobalBasis;
  typedef typename GlobalBasis::GridView GridView;
  typedef typename GlobalBasis::size_type size_type;
  typedef typename GlobalBasis::MultiIndex MultiIndex;
  typedef typename GridView::template Codim<0>::Entity Element;
  typedef PQ1FunctionSpaceBasisLeafNode<GV> Tree;

  PQ1FunctionSpaceBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(globalBasis)
  {}

  //! Bind to element
  void bind(const Element& e)
  {
    element_ = &e;
    tree_.bind(e);
  }

  const Element& element() const
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  //! Unbind from element - only hint
  void unbind()
  {}

  /**
   * \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  const GlobalBasis* globalBasis_;
  const Element* element_;
  Tree tree_;
};


template<typename GV>
class PQ1FunctionSpaceBasisLeafNode
{
public:
  typedef PQ1FunctionSpaceBasis<GV> GlobalBasis;
  typedef typename GlobalBasis::LocalView LocalView;
  typedef typename GlobalBasis::GridView GridView;
  typedef typename GlobalBasis::size_type size_type;
  typedef typename GlobalBasis::MultiIndex MultiIndex;
  typedef typename GridView::template Codim<0>::Entity Element;

  static const int dim = GV::Grid::dimension;

protected:
  typedef typename Dune::PQkLocalFiniteElementCache<typename GV::Grid::ctype, double, dim, 1> FiniteElementCache;

public:
  typedef typename FiniteElementCache::FiniteElementType FiniteElement;

  friend LocalView;

  PQ1FunctionSpaceBasisLeafNode(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    finiteElement_(0),
    element_(0)
  {}

//  PQ1FunctionSpaceBasisLeafNode(const LocalView& localView) :
//    localView_(&localView)
//  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

  //! size of subtree rooted in this node (element-local)
  size_type subTreeSize() const
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_->localBasis().size();
  }

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type maxSubTreeSize() const
  {
    return 1<<dim;
  }

  //! size of complete tree (element-local)
  size_type localSize() const
  {
    // We have localSize==subTreeSize because the tree consist of a single leaf node.
    return subTreeSize();
  }

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const
  {
    return i;
  }

  //! maximum size of complete tree for any element of the global basis
  size_type maxLocalSize() const
  {
    // We have maxLocalSize==maxSubTreeSize because the tree consist of a single leaf node.
    return maxSubTreeSize();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis (pair of multi-indices)
  const MultiIndex globalIndex(size_type i) const
  {
    return { globalBasis_->gridView().indexSet().subIndex(
        *element_,
        finiteElement_->localCoefficients().localKey(i).subEntity(),
        dim) };
  }

  //! Generate multi indices for current subtree into range starting at it
  //! \param it iterator over a container of MultiIndex
  //! \return iterator past the last written element (STL-style)
  template<typename MultiIndexIterator>
  MultiIndexIterator generateMultiIndices(MultiIndexIterator it) const
  {
    size_type size = subTreeSize();
    for(size_type i=0; i<size; ++i)
    {
      (*it) = {globalBasis_->gridView().indexSet().subIndex(
        *element_,
        finiteElement_->localCoefficients().localKey(i).subEntity(),
        dim)};
      ++it;
    }
    return it;
  }

  //! to be discussed
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }


protected:

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
  }

  const GlobalBasis* globalBasis_;
  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1FUNCTIONSPACEBASIS_HH
