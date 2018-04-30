// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BREZZIDOUGLASMARINIBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BREZZIDOUGLASMARINIBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube3d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>

namespace Dune {
namespace Functions {

namespace Impl {

  template<int dim, GeometryType::BasicType basic_type, typename D, typename R, std::size_t k>
  struct BDMLocalInfo
  {
    static_assert((AlwaysFalse<D>::value),"The requested type of BDM element is not implemented, sorry!");
  };

  template<typename D, typename R>
  struct BDMLocalInfo<2,GeometryType::simplex,D,R,1>
  {
    using FiniteElement = BDM1Simplex2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 8;
  };

  template<typename D, typename R>
  struct BDMLocalInfo<2,GeometryType::simplex,D,R,2>
  {
    using FiniteElement = BDM2Simplex2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 8;
  };

  template<typename D, typename R>
  struct BDMLocalInfo<2,GeometryType::cube,D,R,1>
  {
    using FiniteElement = BDM1Cube2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 16;
  };

  template<typename D, typename R>
  struct BDMLocalInfo<2,GeometryType::cube,D,R,2>
  {
    using FiniteElement = BDM2Cube2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 16;
  };

  template<typename D, typename R>
  struct BDMLocalInfo<3,GeometryType::cube,D,R,1>
  {
    using FiniteElement = BDM1Cube3DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 64;
  };

  template<typename GV, int dim, GeometryType::BasicType basic_type, typename D, typename R, std::size_t k>
  class BDMLocalFiniteElementMap
  {
    static const std::size_t Variants = BDMLocalInfo<dim, basic_type, D, R, k>::Variants;

  public:

    using FiniteElement = typename BDMLocalInfo<dim, basic_type, D, R, k>::FiniteElement;

    BDMLocalFiniteElementMap(const GV& gv)
      : gv_(gv), is_(&(gv_.indexSet())), orient_(gv.size(0))
    {
      // create all variants
      for (size_t i = 0; i < Variants; i++)
        variant_[i] = FiniteElement(i);

      // compute orientation for all elements
      // loop once over the grid
      for(const auto& cell : elements(gv))
      {
        unsigned int myId = is_->index(cell);
        orient_[myId] = 0;

        for (const auto& intersection : intersections(gv,cell))
        {
          if (intersection.neighbor() && (is_->index(intersection.outside()) > myId))
            orient_[myId] |= (1 << intersection.indexInInside());
        }
      }
    }

    //! \brief get local basis functions for entity
    template<class EntityType>
    const FiniteElement& find(const EntityType& e) const
    {
      return variant_[orient_[is_->index(e)]];
    }

    private:
      GV gv_;
      std::array<FiniteElement,Variants> variant_;
      const typename GV::IndexSet* is_;
      std::vector<unsigned char> orient_;
  };


} // namespace Impl


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   BrezziDouglasMariniPreBasis
//   BrezziDouglasMariniNodeIndexSet
//   BrezziDouglasMariniNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename ST, typename TP, GeometryType::BasicType basic_type>
class BrezziDouglasMariniNode;

template<typename GV, int k, class MI, class TP, class ST, GeometryType::BasicType basic_type>
class BrezziDouglasMariniNodeIndexSet;

template<typename GV, int k, class MI, class ST, GeometryType::BasicType basic_type>
class BrezziDouglasMariniPreBasis
{
  static const int dim = GV::dimension;
  using FiniteElementMap = typename Impl::BDMLocalFiniteElementMap<GV, dim, basic_type, typename GV::ctype, double, k>;

  template<typename, int, class, class, class, GeometryType::BasicType>
  friend class BrezziDouglasMariniNodeIndexSet;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;

  // Precompute the number of dofs per entity type depending on the entity's codimension
  std::vector<int> dofsPerCodim {dim*(k-1)*3, dim+(k-1)}; // hardcoded

  template<class TP>
  using Node = BrezziDouglasMariniNode<GV, k, size_type, TP, basic_type>;

  template<class TP>
  using IndexSet = BrezziDouglasMariniNodeIndexSet<GV, k, MI, TP, ST, basic_type>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  BrezziDouglasMariniPreBasis(const GridView& gv) :
    gridView_(gv),
    finiteElementMap_(gv)
  { }

  void initializeIndices()
  {
    codimOffset_.resize(3);
    codimOffset_[0] = 0;
    codimOffset_[1] = codimOffset_[0] + dofsPerCodim[0] * gridView_.size(0);
    //if (dim==3) codimOffset_[2] = codimOffset_[1] + dofsPerCodim[1] * gridView_.size(1);
  }

  /* \brief Update the stored grid view, to be called if the grid has changed */
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp, &finiteElementMap_};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    return dofsPerCodim[0] * gridView_.size(0) + dofsPerCodim[1] * gridView_.size(1); // only 2d
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    // The implementation currently only supports grids with a single element type.
    // We can therefore return the actual number of dofs here.
    GeometryType elementType = *(gridView_.indexSet().types(0).begin());
    size_t numFaces = ReferenceElements<double,dim>::general(elementType).size(1);
    return dofsPerCodim[0] + dofsPerCodim[1] * numFaces;
  }

protected:
  const GridView gridView_;
  std::vector<size_t> codimOffset_;
  FiniteElementMap finiteElementMap_;
};



template<typename GV, int k, typename ST, typename TP, GeometryType::BasicType basic_type>
class BrezziDouglasMariniNode :
  public LeafBasisNode<ST, TP>
{
  static const int dim = GV::dimension;

  using Base = LeafBasisNode<ST,TP>;

public:

  using size_type = ST;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElementMap = typename Impl::BDMLocalFiniteElementMap<GV, dim, basic_type, typename GV::ctype, double, k>;
  using FiniteElement = typename FiniteElementMap::FiniteElement;

  BrezziDouglasMariniNode(const TreePath& treePath, const FiniteElementMap* finiteElementMap) :
    Base(treePath),
    finiteElement_(nullptr),
    element_(nullptr),
    finiteElementMap_(finiteElementMap)
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
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(finiteElementMap_->find(*element_));
    this->setSize(finiteElement_->size());
  }

protected:

  const FiniteElement* finiteElement_;
  const Element* element_;
  const FiniteElementMap* finiteElementMap_;
};



template<typename GV, int k, class MI, class TP, class ST, GeometryType::BasicType basic_type>
class BrezziDouglasMariniNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = BrezziDouglasMariniPreBasis<GV, k, MI, ST, basic_type>;

  using Node = typename PreBasis::template Node<TP>;

  BrezziDouglasMariniNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    node_ = &node;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return node_->finiteElement().size();
  }
#if 0
  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  //! Assume dim \in \lbrace 2, 3 \rbrace.
  MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    // The dimension of the entity that the current dof is related to
    size_t subentity = localKey.subEntity();
    size_t codim = localKey.codim();

    // Throw if Element is no cube or simplex
    if (not(basic_type==GeometryType::BasicType::cube and element.type().isCube()) and
        not(basic_type==GeometryType::BasicType::simplex and element.type().isSimplex())) DUNE_THROW(Dune::NotImplemented, "BDMNodalBasis only implemented for cube and simplex elements.");

    if (not(codim==0 or codim==1))
      DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the BDMThomasBasis");

    return { preBasis_->CodimOffset_[codim] +
             preBasis_->dofsPerCodim[codim] * gridIndexSet.subIndex(element, subentity, codim) + localKey.index() };
  }
#endif

  /**
   * \brief Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
   *
   * This assume dim \in \lbrace 2, 3 \rbrace.
   */
  template<typename It>
  It indices(It it) const
  {
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    // throw if element is not of predefined type
    if (not(basic_type==GeometryType::BasicType::cube and element.type().isCube()) and
        not(basic_type==GeometryType::BasicType::simplex and element.type().isSimplex())) DUNE_THROW(Dune::NotImplemented, "BrezziDouglasMariniBasis only implemented for cube and simplex elements.");

    for(std::size_t i=0, end=size(); i<end; ++i, ++it)
    {
      Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);

      // The dimension of the entity that the current dof is related to
      size_t subentity = localKey.subEntity();
      size_t codim = localKey.codim();

      *it = { preBasis_->codimOffset_[codim] +
             preBasis_->dofsPerCodim[codim] * gridIndexSet.subIndex(element, subentity, codim) + localKey.index() };
    }

    return it;
  }

protected:
  const PreBasis* preBasis_;
  const Node* node_;
};



namespace BasisFactory {

namespace Imp {

template<std::size_t k, GeometryType::BasicType basic_type>
class BrezziDouglasMariniPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView, class size_type=std::size_t>
  auto makePreBasis(const GridView& gridView) const
    -> BrezziDouglasMariniPreBasis<GridView, k, MultiIndex, size_type, basic_type>
  {
    return {gridView};
  }
};

} // end namespace BasisFactory::Imp

template<std::size_t k, GeometryType::BasicType basic_type>
Imp::BrezziDouglasMariniPreBasisFactory<k, basic_type> brezziDouglasMarini()
{
  return{};
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order BDM finite element space on simplex and cube grids
 *
 * TODO
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k, GeometryType::BasicType basic_type, class ST = std::size_t>
using BrezziDouglasMariniBasis = DefaultGlobalBasis<BrezziDouglasMariniPreBasis<GV, k, FlatMultiIndex<ST>, ST, basic_type> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BREZZIDOUGLASMARINIBASIS_HH
