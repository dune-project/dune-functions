// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDLAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDLAGRANGEBASIS_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

#include <dune/localfunctions/refined.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>


namespace Dune {
namespace Functions {

template<typename GV, int k, typename R>
class RefinedLagrangeNode;

/**
 * \brief A pre-basis for a refined Lagrange bases
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions
 * \tparam R   Range field-type used for shape function values
 *
 * \note This only works for simplex grids.
 */
template <typename GV, int k, typename R = double>
class RefinedLagrangePreBasis :
  public LeafPreBasisMapperMixIn< GV >
{
  using Base = LeafPreBasisMapperMixIn< GV >;

  static const int dim = GV::dimension;

  // refined basis only implemented for P0 and P1
  static_assert(k == 0 || k == 1);

  // the layout is defined in terms of a MCMGLayout specialized for k == 0 or 1
  static MCMGLayout dofLayout()
  {
    if constexpr(k == 0)
      // a refined P0 basis assigns each element 2^dim DOFs
      return [](GeometryType gt, int) -> size_t {
        return (gt.dim() == dim) ? (1 << dim) : 0;
      };
    else if constexpr(k == 1)
      // a refined P1 basis has the same layout as a P2 basis
      return [](GeometryType gt, int) -> size_t {
        return Dune::binomial(int(k),int(gt.dim()));
      };
    else
      DUNE_THROW(Dune::NotImplemented,
        "Refined basis not implemented for higher-order Lagrange (k>=2) elements.");
  }

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type of the refined Lagrange tree node
  using Node = RefinedLagrangeNode<GV, k, R>;

  /**
   * \brief Constructor for a given grid view object.
   *
   * \param gv The GridView the basis is defined on.
   * \throws Dune::NotImplemented If an element of type !simplex is found.
   */
  RefinedLagrangePreBasis (const GridView& gv)
    : Base(gv, dofLayout())
  {
    for (auto gt : gv.indexSet().types(0)) {
      if (!gt.isSimplex())
        DUNE_THROW(Dune::NotImplemented,
          "Refined Lagrange basis only implemented for simplex grids.");
    }
  }

  //! Create tree node
  Node makeNode () const
  {
    return Node{};
  }

  /**
   * \brief Polynomial order used in the local Lagrange finite-elements.
   *
   * \note The local function is of order `k` only in subdomains of the element.
   *       It might be necessary to use a subdivided quadrature rule for
   *       integration.
   */
  static constexpr unsigned int order()
  {
    return k;
  }
};



template <typename GV, int k, typename R>
class RefinedLagrangeNode
  : public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

  // refined basis only implemented for P0 and P1
  static_assert(k == 0 || k == 1);

public:
  //! Type of the element in the GridView
  using Element = typename GV::template Codim<0>::Entity;

  //! Type of the local finite-element
  using FiniteElement = std::conditional_t<(k==0),
    Dune::RefinedP0LocalFiniteElement<typename GV::ctype,R,dim>,
    Dune::RefinedP1LocalFiniteElement<typename GV::ctype,R,dim>>;

  /**
   * \brief The default constructor initializes all members to their default.
   *
   * The constructor default constructs the local finite-element and sets the
   * element pointer to `nullptr`, meaning that the node is not bound to any
   * element yet.
   *
   * \note Before the node can be used it needs to be bound to an element.
   **/
  RefinedLagrangeNode ()
    : finiteElement_{}
    , element_(nullptr)
  {}

  /**
   * \brief Return current element.
   * The behavior is undefined if the node is not bound to any element.
   */
  const Element& element () const
  {
    return *element_;
  }

  /**
   * \brief Return the LocalFiniteElement for the element we are bound to.
   *
   * The LocalFiniteElement implements the corresponding interfaces of the
   * dune-localfunctions module.
   */
  const FiniteElement& finiteElement () const
  {
    return finiteElement_;
  }

  //! Bind the node to the element `e`.
  void bind (const Element& e)
  {
    element_ = &e;
    this->setSize(finiteElement_.size());
  }

  /**
   * \brief Polynomial order used in the local Lagrange finite-elements in
   * subdomains of the element.
   */
  static constexpr unsigned int order()
  {
    return k;
  }

protected:
  const FiniteElement finiteElement_;
  const Element* element_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a  RefinedLagrange pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 */
template <int k, typename R=double>
auto refinedLagrange ()
{
  return [](const auto& gridView) {
    return RefinedLagrangePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
}

} // end namespace BasisFactory


/** \brief Nodal basis of a continuous Lagrange finite-element space on a uniformly refined simplex element
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for simplex grids.
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of RefinedLagrangeBasis.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 * \tparam R The range type of the local basis
 */
template <typename GV, int k, typename R=double>
using RefinedLagrangeBasis = DefaultGlobalBasis<RefinedLagrangePreBasis<GV,k,R> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDLAGRANGEBASIS_HH
