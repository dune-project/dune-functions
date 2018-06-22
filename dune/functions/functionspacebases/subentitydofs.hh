// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBENTITYDOFS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBENTITYDOFS_HH

#include <vector>
#include <array>
#include <bitset>

#include <dune/geometry/referenceelements.hh>
#include <dune/typetree/traversal.hh>



namespace Dune {
namespace Functions {



/**
 * \brief Range of DOFs associated to sub-entity
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This provides a range of DOFs associated to a given sub-entity
 * by its LocalKeys. In order to use this, it has to be bound to
 * a bound LocalView and a sub-entity. The latter is either encoded
 * by its local index and codimension wrt the element the LocalView
 * is bound to, or by an intersection having this element as inside.
 *
 * After beeing bound this class can be used as range of local indices
 * of those DOFs that are associated to the sub-entity or sub-sub-entities
 * of this sub-entity. Furthermore it allows to ask if a given local index
 * is contained in this range.
 *
 * Notice that the class itself does only depend on the GridView
 * because it needs to allocate some dimension-dependent containers
 * for caching the information computed during bind.
 *
 * \tparam GridView The GridView LocalViews should act on
 */
template<class GridView>
class SubEntityDOFs
{
  static const int dim = GridView::dimension;

  // Caution: This may be insufficient for higher dimensions.
  // However, its enough for storing all standard geometry types
  // up to dim=4 because largest sub-entity count is the number
  // of edges in the tessaract which is 32 edges. Furthermore
  // std::bitset<k> will have size 32 for k<=32 anyway.
  // On 64 bit systems we may also use 64 for the same reason.
  static const int maxSubEntityCount = 32;

  // Fill a table to mark all sub-sub-entities contained in the given
  // sub-entity. The table has one row for each codimension and
  // one column for each sub-entity of this codimension. Rows are
  // indexed with dim-codim.
  //
  // There are some cases where checking 'is contained' is trivial
  // which are thus not stored in the table:
  // * For subEntityCodim=0 that's true for all sub-entities.
  // * For subEntityCodim>codim that's false for all sub-entities.
  //
  // Hence we don't store values for subEntityCodim=0 at all
  // and otherwise only values codim>subEntityCodim.
  //
  void fillIsInSubEntityTable(const Dune::GeometryType& type, std::size_t subEntityIndex, std::size_t subEntityCodim)
  {
    using Dune::referenceElement;
    auto re = referenceElement(double(), type, Dune::Dim<dim>());

    for(auto&& row : isInSubEntity_)
      row.reset();

    // We'd like to loop over all sub-entities, but this is not
    // supported by the reference elements yet:
    // forEachSubEntity(re, subEntityIndex, subEntityCodim, [&](auto i, auto codim) {
    //   isInSubEntity_[dim-codim][i] = true;
    // });

    // For codim<subEntityCodim no sub-entity is contained and
    // for codim=subEntityCodim only the subEntityIndex-th sub-entity
    // is contained. Hence we treat the latter case explicitly and
    // only loop over codim>=subEntityCodim+1.
    isInSubEntity_[dim-subEntityCodim][subEntityIndex] = true;
    for(std::size_t codim=subEntityCodim+1; codim<=dim; ++codim)
      for(std::size_t i=0; i<static_cast<std::size_t>(re.size(subEntityIndex, subEntityCodim, codim)); ++i)
        isInSubEntity_[dim-codim][re.subEntity(subEntityIndex, subEntityCodim, i, codim)] = true;
  }

public:

  /**
   * \brief Bind SubEntityDOFs object to LocalView and sub-entity
   *
   * Notice that this method will pre-compute and cache the contained
   * local indices as well as a look-up table for checking if a local
   * index is contained.
   * The sub-entity is encoded as (index,codim) with respect to
   * the element the LocalView is bound to.
   * In order to be able to bind and use the range in a single
   * expression this returns *this for convenience.
   *
   * \param localView A bound LocalView to bind to
   * \param subEntityIndex Index of sub-entity in localView.element()
   * \param subEntityCodim Codimension of sub-entity in localView.element()
   * \returns *this for convenience
   */
  template<class LocalView>
  SubEntityDOFs& bind(const LocalView& localView, std::size_t subEntityIndex, std::size_t subEntityCodim)
  {
    const auto& tree = localView.tree();

    // fill vector with local indices of all DOFs contained in subentity
    containedDOFs_.clear();
    dofIsContained_.clear();
    dofIsContained_.resize(localView.size(), false);
    if (subEntityCodim == 0)
      for (std::size_t i=0; i<tree.size(); ++i)
      {
        containedDOFs_.push_back(tree.localIndex(i));
        dofIsContained_[tree.localIndex(i)] = true;
      }
    else
    {
      fillIsInSubEntityTable(localView.element().type(), subEntityIndex, subEntityCodim);

      Dune::TypeTree::forEachLeafNode(tree, [&, dim=dim](auto&& node, auto&& treePath) {
        const auto& localCoefficients = node.finiteElement().localCoefficients();
        std::size_t localSize = localCoefficients.size();
        for(std::size_t i=0; i<localSize; ++i)
        {
          auto localKey = localCoefficients.localKey(i);
          if ((localKey.codim() >= subEntityCodim) and (isInSubEntity_[dim-localKey.codim()][localKey.subEntity()]))
          {
            containedDOFs_.push_back(node.localIndex(i));
            dofIsContained_[tree.localIndex(i)] = true;
          }
        }
      });
    }
    return *this;
  }

  /**
   * \brief Bind SubEntityDOFs object to LocalView and sub-entity
   *
   * Notice that this method will pre-compute and cache the contained
   * local indices as well as a look-up table for checking if a local
   * index is contained.
   * The sub-entity is encoded as intersection having localView.element()
   * as inside. This is equivalent to bind(localView, intersection.indexInInside(),1),
   * In order to be able to bind and use the range in a single
   * expression this returns *this for convenience.
   *
   * \param localView A bound LocalView to bind to
   * \param intersection An Intersection encoding the sub-entity to bind to
   * \returns *this for convenience
   */
  template<class LocalView, class Intersection>
  SubEntityDOFs& bind(const LocalView& localView, const Intersection& intersection)
  {
    return bind(localView, intersection.indexInInside(), 1);
  }

  //! Create begin iterator for access to range of contained local indices
  auto begin() const
  {
    return containedDOFs_.cbegin();
  }

  //! Create end iterator for access to range of contained local indices
  auto end() const
  {
    return containedDOFs_.cend();
  }

  //! Return number of contained DOFs
  auto size() const
  {
    return containedDOFs_.size();
  }

  //! Return i-th entry of the range of contained local indices
  decltype(auto) operator[](std::size_t i) const
  {
    return containedDOFs_[i];
  }

  //! Check if given local index is contained in this range of DOFs
  bool contains(std::size_t localIndex) const
  {
    return dofIsContained_[localIndex];
  }

private:

  std::array< std::bitset<maxSubEntityCount>, dim> isInSubEntity_;
  std::vector<std::size_t> containedDOFs_;
  std::vector<bool> dofIsContained_;
};



/**
 * \brief Create SubEntityDOFs object
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * The only requirement on the passed argument t
 * is, that its type T has to provide a typedef
 * T::GridView. Hence t it can be a GlobalBasis or
 * a LocalView.
 *
 * \param t A GlobalBasis or a LocalView
 */
template<class T>
auto subEntityDOFs(const T& t)
{
  return SubEntityDOFs<typename T::GridView>{};
}



/**
 * \brief Create bound SubEntityDOFs object
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This creates a SubEntityDOFs object and binds
 * it to the given LocalView and sub-entity.
 *
 * Notice that the SubEntityDOFs object will allocate
 * some internal buffers. For efficiency reasons you
 * should thus prefer to first create a SubEntityDOFs
 * object and then bind it to each element you want to
 * process instead of creating a new bound SubEntityDOFs
 * object for each element.
 *
 * \param localView A bound LocalView to bind to
 * \param subEntityIndex Index of sub-entity in localView.element()
 * \param subEntityCodim Codimension of sub-entity in localView.element()
 */
template<class LocalView>
auto subEntityDOFs(const LocalView& localView, std::size_t subEntityIndex, std::size_t subEntityCodim)
{
  using GridView = typename LocalView::GridView;
  SubEntityDOFs<GridView> subEntityDOFs;
  subEntityDOFs.bind(localView, subEntityIndex, subEntityCodim);
  return subEntityDOFs;
}



/**
 * \brief Create bound SubEntityDOFs object
 *
 * \ingroup FunctionSpaceBasesUtilities
 *
 * This creates a SubEntityDOFs object and binds
 * it to the given LocalView and intersection.
 *
 * Notice that the SubEntityDOFs object will allocate
 * some internal buffers. For efficiency reasons you
 * should thus prefer to first create a SubEntityDOFs
 * object and then bind it to each element you want to
 * process instead of creating a new bound SubEntityDOFs
 * object for each element.
 *
 * \param localView A bound LocalView to bind to
 * \param intersection An Intersection encoding the sub-entity to bind to
 */
template<class LocalView, class Intersection>
auto subEntityDOFs(const LocalView& localView, const Intersection& intersection)
{
  using GridView = typename LocalView::GridView;
  SubEntityDOFs<GridView> subEntityDOFs;
  subEntityDOFs.bind(localView, intersection);
  return subEntityDOFs;
}



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SUBENTITYDOFS_HH
