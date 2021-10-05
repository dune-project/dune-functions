// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BREZZIDOUGLASMARINIBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BREZZIDOUGLASMARINIBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube3d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh>

#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

namespace Dune {
namespace Functions {

namespace Impl {

  template<int dim, typename D, typename R, std::size_t k>
  struct BDMSimplexLocalInfo
  {
    static_assert((AlwaysFalse<D>::value),"The requested type of BDM element is not implemented, sorry!");
  };

  template<typename D, typename R>
  struct BDMSimplexLocalInfo<2,D,R,1>
  {
    using FiniteElement = BDM1Simplex2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 8;
  };

  template<typename D, typename R>
  struct BDMSimplexLocalInfo<2,D,R,2>
  {
    using FiniteElement = BDM2Simplex2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 8;
  };

  template<int dim, typename D, typename R, std::size_t k>
  struct BDMCubeLocalInfo
  {
    static_assert((AlwaysFalse<D>::value),"The requested type of BDM element is not implemented, sorry!");
  };

  template<typename D, typename R>
  struct BDMCubeLocalInfo<2,D,R,1>
  {
    using FiniteElement = BDM1Cube2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 16;
  };

  template<typename D, typename R>
  struct BDMCubeLocalInfo<2,D,R,2>
  {
    using FiniteElement = BDM2Cube2DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 16;
  };

  template<typename D, typename R>
  struct BDMCubeLocalInfo<3,D,R,1>
  {
    using FiniteElement = BDM1Cube3DLocalFiniteElement<D,R>;
    static const std::size_t Variants = 64;
  };

  template<typename GV, int dim, typename R, std::size_t k>
  class BDMLocalFiniteElementMap
  {
    using D = typename GV::ctype;
    using CubeFiniteElement    = typename BDMCubeLocalInfo<dim, D, R, k>::FiniteElement;
    using SimplexFiniteElement = typename BDMSimplexLocalInfo<dim, D, R, k>::FiniteElement;

  public:

    using T = LocalBasisTraits<D, dim, FieldVector<D,dim>, R, dim, FieldVector<R,dim>, FieldMatrix<D,dim,dim> >;
    using FiniteElement = LocalFiniteElementVirtualInterface<T>;

    BDMLocalFiniteElementMap(const GV& gv)
      : is_(&(gv.indexSet())), orient_(gv.size(0))
    {
      cubeVariant_.resize(BDMCubeLocalInfo<dim, D, R, k>::Variants);
      simplexVariant_.resize(BDMSimplexLocalInfo<dim, D, R, k>::Variants);

      // create all variants
      for (size_t i = 0; i < cubeVariant_.size(); i++)
        cubeVariant_[i] = std::make_unique<LocalFiniteElementVirtualImp<CubeFiniteElement> >(CubeFiniteElement(i));

      for (size_t i = 0; i < simplexVariant_.size(); i++)
        simplexVariant_[i] = std::make_unique<LocalFiniteElementVirtualImp<SimplexFiniteElement> >(SimplexFiniteElement(i));

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
      if (e.type().isCube())
        return *cubeVariant_[orient_[is_->index(e)]];
      else
        return *simplexVariant_[orient_[is_->index(e)]];
    }

    private:
      std::vector<std::unique_ptr<LocalFiniteElementVirtualImp<CubeFiniteElement> > > cubeVariant_;
      std::vector<std::unique_ptr<LocalFiniteElementVirtualImp<SimplexFiniteElement> > > simplexVariant_;
      const typename GV::IndexSet* is_;
      std::vector<unsigned char> orient_;
  };


} // namespace Impl


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   BrezziDouglasMariniPreBasis
//   BrezziDouglasMariniNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k>
class BrezziDouglasMariniNode;

template<typename GV, int k>
class BrezziDouglasMariniPreBasis
{
  static const int dim = GV::dimension;
  using FiniteElementMap = typename Impl::BDMLocalFiniteElementMap<GV, dim, double, k>;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  using Node = BrezziDouglasMariniNode<GV, k>;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  /** \brief Constructor for a given grid view object */
  BrezziDouglasMariniPreBasis(const GridView& gv) :
    gridView_(gv),
    finiteElementMap_(gv)
  {
    // There is no inherent reason why the basis shouldn't work for grids with more than one
    // element types.  Somebody simply has to sit down and implement the missing bits.
    if (gv.indexSet().types(0).size() > 1)
      DUNE_THROW(Dune::NotImplemented, "Brezzi-Douglas-Marini basis is only implemented for grids with a single element type");
  }

  void initializeIndices()
  {
    codimOffset_[0] = 0;
    codimOffset_[1] = codimOffset_[0] + dofsPerCodim_[0] * gridView_.size(0);
    //if (dim==3) codimOffset_[2] = codimOffset_[1] + dofsPerCodim[1] * gridView_.size(1);
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  /* \brief Update the stored grid view, to be called if the grid has changed */
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{&finiteElementMap_};
  }

  size_type size() const
  {
    return dofsPerCodim_[0] * gridView_.size(0) + dofsPerCodim_[1] * gridView_.size(1); // only 2d
  }

  //! Return number possible values for next position in multi index
  template<class SizePrefix>
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
    return dofsPerCodim_[0] + dofsPerCodim_[1] * numFaces;
  }

  /**
   * \brief Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
   *
   * This assume dim \in \lbrace 2, 3 \rbrace.
   */
  template<typename It>
  It indices(const Node& node, It it) const
  {
    const auto& gridIndexSet = gridView().indexSet();
    const auto& element = node.element();

    // throw if element is not of predefined type
    if (not(element.type().isCube()) and not(element.type().isSimplex()))
      DUNE_THROW(Dune::NotImplemented, "BrezziDouglasMariniBasis only implemented for cube and simplex elements.");

    for(std::size_t i=0, end=node.size(); i<end; ++i, ++it)
    {
      Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);

      // The dimension of the entity that the current dof is related to
      size_t subentity = localKey.subEntity();
      size_t codim = localKey.codim();

      *it = { codimOffset_[codim] +
             dofsPerCodim_[codim] * gridIndexSet.subIndex(element, subentity, codim) + localKey.index() };
    }

    return it;
  }

protected:
  GridView gridView_;
  std::array<size_t,dim+1> codimOffset_;
  FiniteElementMap finiteElementMap_;
  // Number of dofs per entity type depending on the entity's codimension and type
  std::array<int,2> dofsPerCodim_ {{dim*(k-1)*3, dim+(k-1)}};
};



template<typename GV, int k>
class BrezziDouglasMariniNode :
  public LeafBasisNode
{
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElementMap = typename Impl::BDMLocalFiniteElementMap<GV, dim, double, k>;
  using FiniteElement = Impl::GlobalValuedLocalFiniteElement<Impl::ContravariantPiolaTransformator,
                                                             typename FiniteElementMap::FiniteElement,
                                                             Element>;

  BrezziDouglasMariniNode(const FiniteElementMap* finiteElementMap) :
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
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_.bind((finiteElementMap_->find(*element_)), e);
    this->setSize(finiteElement_.size());
  }

protected:

  FiniteElement finiteElement_;
  const Element* element_;
  const FiniteElementMap* finiteElementMap_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a Brezzi-Douglas-Marini pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k Order of the Brezzi-Douglas-Marini element
 */
template<std::size_t k>
auto brezziDouglasMarini()
{
  return [](const auto& gridView) {
    return BrezziDouglasMariniPreBasis<std::decay_t<decltype(gridView)>, k>(gridView);
  };
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order BDM finite element space on simplex and cube grids
 *
 * TODO: Fix this for grids with more than one element type
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using BrezziDouglasMariniBasis = DefaultGlobalBasis<BrezziDouglasMariniPreBasis<GV, k> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BREZZIDOUGLASMARINIBASIS_HH
