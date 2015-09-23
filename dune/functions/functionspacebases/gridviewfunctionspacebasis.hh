// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH

#include <dune/typetree/leafnode.hh>
#include <dune/functions/functionspacebases/nodes.hh>

namespace Dune {
namespace Functions {



template<typename E, typename FE, typename ST, typename TP>
class GridFunctionSpaceBasisLeafNodeInterface :
  public LeafBasisNode<ST, TP>
{
  using Base = LeafBasisNode<ST, TP>;
public:
  typedef ST size_type;
  typedef E Element;
  typedef FE FiniteElement;

  using TreePath = typename Base::TreePath;

  GridFunctionSpaceBasisLeafNodeInterface(TreePath treePath = TreePath()) :
    Base(treePath)
  {}


  //! Return current element, throw if unbound
  virtual const Element& element() const = 0;

  virtual const FiniteElement& finiteElement() const = 0;

  //! size of subtree rooted in this node (element-local)
  virtual size_type size() const = 0;

};

template<typename GV, typename LV, typename IS, typename MI>
class GridViewFunctionSpaceBasis
{
public:

  typedef GV GridView;
  typedef std::size_t size_type;
  typedef LV LocalView;
  typedef IS IndexSet;
  typedef MI MultiIndex;

  /** \brief Obtain the grid view that the basis is defined on
   */
  virtual const GridView& gridView() const = 0;

  /**
   * \brief Return local view for basis
   *
   * Perhaps we must move the construction outside
   * of the global basis in order to calm the compiler
   * when instantiating the TMP constructing the local view.
   */
  virtual LocalView localView() const = 0;

};

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH
