// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ALGORITHMS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ALGORITHMS_HH



#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>



namespace Dune {
namespace Functions {



namespace TreeAlgorithms {



namespace Impl {

  // A TreeVisitor using plain callbacks to implement the
  // visitor pattern.
  template<class PreFunc, class LeafFunc, class PostFunc>
  struct CallbackVisitor :
    public Dune::TypeTree::TreeVisitor,
    public Dune::TypeTree::DynamicTraversal
  {
    public:
    CallbackVisitor(PreFunc& preFunc, LeafFunc& leafFunc, PostFunc& postFunc) :
      preFunc_(preFunc),
      leafFunc_(leafFunc),
      postFunc_(postFunc)
    {}

    template<typename Node, typename TreePath>
    void pre(Node&& node, TreePath treePath)
    {
      preFunc_(node, treePath);
    }

    template<typename Node, typename TreePath>
    void leaf(Node&& node, TreePath treePath)
    {
      leafFunc_(node, treePath);
    }

    template<typename Node, typename TreePath>
    void post(Node&& node, TreePath treePath)
    {
      postFunc_(node, treePath);
    }

  private:
    PreFunc& preFunc_;
    LeafFunc& leafFunc_;
    PostFunc& postFunc_;
  };

  // Convinience function for creating a CallbackVisitor
  template<class PreFunc, class LeafFunc, class PostFunc>
  auto callbackVisitor(PreFunc& preFunc, LeafFunc& leafFunc, PostFunc& postFunc)
  {
    return CallbackVisitor<PreFunc, LeafFunc, PostFunc>(preFunc, leafFunc, postFunc);
  }

} // namespace Impl



/**
 * \brief Traverse tree and visit each node
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * All passed callback functions are called with the
 * node and corresponding treepath as arguments.
 *
 * \param tree The tree to traverse
 * \param preFunc This function is called for each inner node before visiting its children
 * \param leafFunc This function is called for each leaf node
 * \param postFunc This function is called for each inner node after visiting its children
 */
template<class Tree, class PreFunc, class LeafFunc, class PostFunc>
auto forEachNode(Tree&& tree, PreFunc&& preFunc, LeafFunc&& leafFunc, PostFunc&& postFunc)
{
  Dune::TypeTree::applyToTree(tree, Impl::callbackVisitor(preFunc, leafFunc, postFunc));
}

/**
 * \brief Traverse tree and visit each node
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * All passed callback functions are called with the
 * node and corresponding treepath as arguments.
 *
 * \param tree The tree to traverse
 * \param innerFunc This function is called for each inner node before visiting its children
 * \param leafFunc This function is called for each leaf node
 */
template<class Tree, class InnerFunc, class LeafFunc>
auto forEachNode(Tree&& tree, InnerFunc&& innerFunc, LeafFunc&& leafFunc)
{
  auto nop = [](auto&&... args) {};
  forEachNode(tree, innerFunc, leafFunc, nop);
}

/**
 * \brief Traverse tree and visit each node
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * The passed callback function is called with the
 * node and corresponding treepath as arguments.
 *
 * \param tree The tree to traverse
 * \param nodeFunc This function is called for each node
 */
template<class Tree, class NodeFunc>
auto forEachNode(Tree&& tree, NodeFunc&& nodeFunc)
{
  forEachNode(tree, nodeFunc, nodeFunc);
}

/**
 * \brief Traverse tree and visit each leaf node
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * The passed callback function is called with the
 * node and corresponding treepath as arguments.
 *
 * \param tree The tree to traverse
 * \param leafFunc This function is called for each leaf node
 */
template<class Tree, class LeafFunc>
auto forEachLeafNode(Tree&& tree, LeafFunc&& leafFunc)
{
  auto nop = [](auto&&... args) {};
  forEachNode(tree, nop, leafFunc, nop);
}



} // end namespace TreeAlgorithms



namespace BasisAlgorithms {



/**
 * \brief Lopp over sub-sub-entities of reference element sub-entity
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * Given a referenceElement any sub-entity can be identified
 * by the pair s=(subEntityIndex,subEntityCodim).  For a given
 * sub-entity s this function loop over s visiting all sub-sub-entities
 * of s with given codimension subSubEntityCodim
 * (wrt the reference element).
 *
 * \param referenceElement The underlying reference element
 * \param subEntityIndex Index of sub-entity to loop over
 * \param subEntityCodim Codimension of sub entity to loop over
 * \param subSubEntityCodim Codimension of sub-sub-entities to vist.
 * \param f A callback which will be called with the index of any visited sub-sub-entity.
 */
template<class ctype, int dim, class F>
void forEachSubEntity(const Dune::ReferenceElement<ctype, dim>& referenceElement, std::size_t subEntityIndex, std::size_t subEntityCodim, std::size_t subSubEntityCodim, F&& f)
{
  for(std::size_t i=0; i<referenceElement.size(subEntityIndex, subEntityCodim, subSubEntityCodim); ++i)
    f(referenceElement.subEntity(subEntityIndex, subEntityCodim, i, subSubEntityCodim));
}



/**
 * \brief Lopp over sub-sub-entities of reference element sub-entity
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * Given a referenceElement any sub-entity can be identified
 * by the pair s=(subEntityIndex,subEntityCodim).  For a given
 * sub-entity s this function loop over s visiting all sub-sub-entities
 * of s. This will consider sub-sub-enitities of any codimension.
 *
 * \param referenceElement The underlying reference element
 * \param subEntityIndex Index of sub-entity to loop over
 * \param subEntityCodim Codimension of sub entity to loop over
 * \param f A callback which will be called with the index and codimension of any visited sub-sub-entity.
 */
template<class ctype, int dim, class F>
void forEachSubEntity(const Dune::ReferenceElement<ctype, dim>& referenceElement, std::size_t subEntityIndex, std::size_t subEntityCodim, F&& f)
{
  for(std::size_t subSubEntityCodim = subEntityCodim; subSubEntityCodim<=dim; ++subSubEntityCodim)
    forEachSubEntity(referenceElement, subEntityIndex, subEntityCodim, subSubEntityCodim, [&](std::size_t i) { f(i, subSubEntityCodim);});
}



/**
 * \brief Loop over elements of grid view for basis
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This loops over the elements in the grid view associated
 * to the basis and provides a local view a local index set
 * for each element.
 *
 * \param basis A function space basis
 * \param f A callback that will be called with element, bound local view, and bound local index set
 */
template<class Basis, class F>
void forEachElement(const Basis& basis, F&& f)
{
  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();
  auto&& gridView = basis.gridView();
  for(auto&& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);
    f(e, localView, localIndexSet);
  }
}



/**
 * \brief Loop over all DOFs in intersection
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This loops over all DOFs of a basis associated to sub-entities
 * of a given intersection.
 *
 * \param intersection The intersection under consideration
 * \param localView A local view of the basis bount to the inside element of the intersection
 * \param f A callback that will be called with the local index of any DOF associated to a sub-entity of the intersection
 */
template<class Intersection, class LocalView, class F>
void forEachIntersectionDOF(const Intersection& intersection, const LocalView& localView, F&& f)
{
  using Element = typename Intersection::Entity;
  static const int dim = Element::dimension;

  auto&& element = intersection.inside();
  const auto& referenceElement = Dune::ReferenceElements<double, dim>::general(element.type());

  auto faceIndex = intersection.indexInInside();

  // Caution: This is currently limited to meshes with dim<=3
  auto isInFace = std::array< std::array<bool,12>, dim-1>();
  for(auto&& row : isInFace)
    row.fill(false);

//  forEachSubEntity(referenceElement, faceIndex, 1, [&](auto i, auto codim) {
//    isInFace[dim-codim][i] = true;
//  });

  for(std::size_t codim=2; codim<=dim; ++codim)
    for(std::size_t i=0; i<referenceElement.size(faceIndex, 1, codim); ++i)
      isInFace[dim-codim][referenceElement.subEntity(faceIndex, 1, i, codim)] = true;

  Dune::TypeTree::forEachLeafNode(localView.tree(), [&](auto&& node, auto&& treePath) {
    const auto& localBasis = node.finiteElement().localBasis();
    const auto& localCoefficients = node.finiteElement().localCoefficients();
    std::size_t localSize = localBasis.size();
    for(std::size_t i=0; i<localSize; ++i)
    {
        auto localKey = localCoefficients.localKey(i);
        if ((localKey.codim() == 1) and (localKey.subEntity() == faceIndex))
          f(node.localIndex(i));
        if ((localKey.codim() > 1) and (isInFace[dim-localKey.codim()][localKey.subEntity()]))
          f(node.localIndex(i));
    }
  });
}



/**
 * \brief Loop over all DOFs on the boundary
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This loops over all DOFs of a basis associated to sub-entities
 * on the boundary.
 *
 * \param basis A function space basis
 * \param f A callback that will be called with local index, local view, local index set, and intersection of the visited boundary DOF
 */
template<class Basis, class F>
void forEachBoundaryDOF(const Basis& basis, F&& f)
{
  forEachElement(basis, [&](auto&& element, auto&& localView, auto&& localIndexSet) {
    for(const auto& intersection: intersections(basis.gridView(), element))
      if (intersection.boundary())
        forEachIntersectionDOF(intersection, localView, [&](std::size_t localIndex) {
            f(localIndex, localView, localIndexSet, intersection);
        });
  });
}



} // namespace BasisAlgorithms
} // namespace Functions
} // namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ALGORITHMS_HH
