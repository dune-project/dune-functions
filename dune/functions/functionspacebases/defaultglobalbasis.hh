// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH

#include <cstddef>
#include <type_traits>
#include <utility>

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/concept.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/defaultlocalview.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>



namespace Dune {
namespace Functions {



/**
 * \brief Global basis for given pre-basis
 *
 * This class implements the interface of a global basis
 * using the details from a given pre-basis. Hence
 * it serves as an example for this interface.
 *
 * If you want to implement your own global basis, it may be
 * better to implement a pre-basis instead. On the one hand
 * this needs less boiler-plate code. On the other hand
 * it makes your implementation composable and thus much
 * more flexible. That is, you can reuse your pre-basis
 * as one part in a larger product space by plugging it
 * e.g. into a CompositePreBasis of PowerPreBasis.
 * The actual global basis for your FooPreBasis is
 * then obtained by using DefaultGlobalBasis<FooPreBasis>.
 *
 * \tparam PB  Pre-basis providing the implementation details
 */
template<class PB>
class DefaultGlobalBasis
{
public:

  //! Pre-basis providing the implementation details
  using PreBasis = PB;

  //! The empty prefix path that identifies the root in the local ansatz tree
  using PrefixPath = TypeTree::HybridTreePath<>;

  //! The grid view that the FE space is defined on
  using GridView = typename PreBasis::GridView;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = typename PreBasis::MultiIndex;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Type of the local view on the restriction of the basis to a single element
  using LocalView = DefaultLocalView<DefaultGlobalBasis<PreBasis>>;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = typename PreBasis::SizePrefix;

  /**
   * \brief Constructor
   *
   * \tparam T Argument list for PreBasis
   * \param t Argument list for PreBasis
   *
   * This will forward all arguments to the constructor of PreBasis
   */
  template<class... T,
    disableCopyMove<DefaultGlobalBasis, T...> = 0,
    enableIfConstructible<PreBasis, T...> = 0>
  DefaultGlobalBasis(T&&... t) :
    preBasis_(std::forward<T>(t)...),
    prefixPath_()
  {
    static_assert(models<Concept::PreBasis<GridView>, PreBasis>(), "Type passed to DefaultGlobalBasis does not model the PreBasis concept.");
    preBasis_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return preBasis_.gridView();
  }

  //! Obtain the pre-basis providing the implementation details
  const PreBasis& preBasis() const
  {
    return preBasis_;
  }

  //! Obtain the pre-basis providing the implementation details
  PreBasis& preBasis()
  {
    return preBasis_;
  }

  /**
   * \brief Update the stored grid view
   *
   * This will update the indexing information of the global basis.
   * It must be called if the grid has changed.
   */
  void update(const GridView & gv)
  {
    preBasis_.update(gv);
    preBasis_.initializeIndices();
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return preBasis_.dimension();
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return preBasis_.size();
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return preBasis_.size(prefix);
  }

  //! Return local view for basis
  LocalView localView() const
  {
    return LocalView(*this);
  }

  //! Return *this because we are not embedded in a larger basis
  const DefaultGlobalBasis& rootBasis() const
  {
    return *this;
  }

  //! Return empty path, because this is the root in the local ansatz tree
  const PrefixPath& prefixPath() const
  {
    return prefixPath_;
  }

protected:
  PreBasis preBasis_;
  PrefixPath prefixPath_;
};



namespace BasisFactory {

template<class GridView, class PreBasisFactory>
auto makeBasis(const GridView& gridView, PreBasisFactory&& preBasisFactory)
{
  using RawPreBasisFactory = std::decay_t<PreBasisFactory>;
  using MultiIndex = std::conditional_t<
    (RawPreBasisFactory::requiredMultiIndexSize == 1),
    FlatMultiIndex<std::size_t>,
    Dune::ReservedVector<std::size_t, RawPreBasisFactory::requiredMultiIndexSize>>;
  auto preBasis = preBasisFactory.template makePreBasis<MultiIndex>(gridView);
  using PreBasis = std::decay_t<decltype(preBasis)>;

  return DefaultGlobalBasis<PreBasis>(std::move(preBasis));
}

template<class MultiIndex, class GridView, class PreBasisFactory>
auto makeBasis(const GridView& gridView, PreBasisFactory&& preBasisFactory)
{
  auto preBasis = preBasisFactory.template makePreBasis<MultiIndex>(gridView);
  using PreBasis = std::decay_t<decltype(preBasis)>;

  return DefaultGlobalBasis<PreBasis>(std::move(preBasis));
}

} // end namespace BasisFactory

// Backward compatibility
namespace BasisBuilder {

  using namespace BasisFactory;

}

namespace Impl {

  /** \brief Transform values of affine families of finite elements.
   *
   * For affine families, the required transformation is the identity.
   *
   * \param valuesLocal The shape function values to transform
   * \param local The position in the reference element where the shape functions have been evaluated
   * \param geometry The grid elements in world coordinates where the values should be transformed to
   */
  struct EmptyTransformator
  {
    template<typename Values, typename LocalCoordinate, typename Geometry>
    auto apply(Values&& valuesLocal,
      const LocalCoordinate& xi,
      const Geometry& geometry)
    {
      return std::move(valuesLocal);
    }
  };

  /** \brief Transforms shape function values and derivatives from reference element coordinates
   *   to world coordinates using the Piola transform
   */
  struct PiolaTransformator
  {
    template<typename Values, typename LocalCoordinate, typename Geometry>
    auto apply(Values&& valuesLocal,
      const LocalCoordinate& xi,
      const Geometry& geometry)
    {
      auto jacobianTransposed = geometry.jacobianTransposed(xi);
      auto integrationElement = geometry.integrationElement(xi);

      for (auto& value : valuesLocal)
      {
        auto tmp = value;
        jacobianTransposed.mtv(tmp, value);
        value /= integrationElement;
      }

      return std::move(valuesLocal);
    }
  };

/** \brief A class that transforms shape function values and derivatives
 *   from reference element coordinates to world coordinates
 *
 * This default implementation contains the case of affine families of
 * finite elements, where the shape function values do not need any transformation
 * at all, and the derivates need to be multiplied from the left with the
 * transposed inverse Jacobian.
 *
 * Function space basis implementations that require a different transformation
 * need to specialize this class.
 */
template <class BasisTreeNode>
auto getToGlobalTransformator(BasisTreeNode)
{
  return EmptyTransformator();
}

} // end namespace Impl


} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTGLOBALBASIS_HH
