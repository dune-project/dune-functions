// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DYNAMIC_QK_NODE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DYNAMIC_QK_NODE_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/dynamicdgbasis.hh>

namespace Dune {
namespace Functions {


/* Example implementation of a dynamic DG node. This assumes that the geometries are the same
 * for all elements. If one does not have this property (or have other assumptions on your nodes)
 * one can implement another Node based on this template */
template<typename GV, typename TP, typename Cache>
class DynamicOrderNode :
  public LeafBasisNode<std::size_t, TP>
{
  static const int dim = GV::dimension;

  using Base = LeafBasisNode<std::size_t,TP>;
  using FiniteElementCache = Cache;
  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
  using DegreeMap = std::vector<int>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  DynamicOrderNode(const TreePath& treePath, const Mapper* m, const DegreeMap* dgm, const std::vector<size_type>* offSets) :
    Base(treePath),
    finiteElement_(nullptr),
    element_(nullptr),
    mcmgMapper_(m),
    degreeMap_(dgm),
    offSets_(offSets)
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
    return cache_.get(degree_);
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    auto i = mcmgMapper_->index(e);

    degree_ = degreeMap_->at(i);
    offSet_ = offSets_->at(i);

    finiteElement_ = &(cache_.get(degree_));
    this->setSize(finiteElement_->size());
  }

  size_type offSet() const {
    return offSet_;
  }

protected:

  /* friend-ify all the things! Maybe there's a better way... TODO */
  template<typename, typename>
  friend class DynamicDGPreBasis;

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
  const Mapper* mcmgMapper_;
  const DegreeMap* degreeMap_;
  int degree_ = 1;
  const std::vector<size_type>* offSets_;
  size_type offSet_;
};

} // end namespace Functions
} // end namespace Dune

#endif
