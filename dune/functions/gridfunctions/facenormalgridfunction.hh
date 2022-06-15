// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_FACENORMALGRIDFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_FACENORMALGRIDFUNCTION_HH

#include <type_traits>
#include <optional>

#include <dune/common/exceptions.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>


namespace Dune::Functions {

namespace Impl {

// Compute closest face to point
template<class ReferenceElement, class Coordinate>
auto closestFaceIndex(const ReferenceElement& re, const Coordinate& x)
{
  auto closestFaceIndex = decltype(re.subEntity(0,1,0,1)){};
  double closestFaceDistance = std::numeric_limits<double>::max();
  for(auto&& faceIndex : Dune::range(re.size(1)))
  {
    // For a face unit outer normal consider the orthogonal projection
    // Px = x + <c-x,n>*n into the face. Then the distance to the face
    // is given by |x-Px| = |<c-x,n>||n| = <c-x,n>.
    auto normal = re.integrationOuterNormal(faceIndex);
    normal /= normal.two_norm();
    auto c = re.position(faceIndex,1);
    c -= x;
    auto faceDistance = (c*normal);
    if (faceDistance<closestFaceDistance)
    {
      closestFaceDistance = faceDistance;
      closestFaceIndex = faceIndex;
    }
  }
  return closestFaceIndex;
}

} // end namespace Impl




/**
 * \brief Grid function implementing the piecewise element face normal.
 *
 * \ingroup FunctionImplementations
 *
 * This function implements the grid interface.
 * When evaluated at a point x inside of an element, it computes
 * the unit outward normal vector of the face closest to this
 * point.
 *
 * \tparam GV  The GridView associated to the grid-function
 */
template<class GV>
class FaceNormalGridFunction
{
public:
  using GridView = GV;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = typename EntitySet::GlobalCoordinate;

private:

  using Traits = Imp::GridFunctionTraits<Range(Domain), EntitySet, DefaultDerivativeTraits, 16>;

  class LocalFunction
  {
    using Geometry = typename Element::Geometry;
    static const int dimension = GV::dimension;
  public:

    /**
     * \brief Bind the local-function to the passed element.
     *
     * This function stores a copy of `element` and its geometry.
     *
     * \b Expects:
     * - The `element` is in the entitySet of the grid-fuction.
     *
     * \b Ensures:
     * - The local-function is bound the the `element`.
     **/
    void bind(const Element& element)
    {
      element_ = element;
      geometry_.emplace(element_.geometry());
    }

    void unbind()
    {
      geometry_.reset();
    }

    /** \brief Return if the local function is bound to a grid element
     */
    bool bound() const
    {
      return static_cast<bool>(geometry_);
    }

    /**
     * \brief Evaluate the local-function in local coordinates `x`.
     *
     * This function computes the unit outward normal vector of the
     * face closest to the point `x` in the bound element.
     *
     * \b Expects:
     * - The local-function is bound to an element.
     **/
    Range operator()(const LocalDomain& x) const
    {
      auto&& re = Dune::referenceElement(*geometry_);
      // Compute reference normal of closest face to given point
      auto face = Impl::closestFaceIndex(re, x);
      auto localNormal = re.integrationOuterNormal(face);

      // Transform reference normal into global unit outer normal using
      // covariant Piola transformation
      auto normal = Range{};
      geometry_->jacobianInverseTransposed(x).mv(localNormal, normal);
      normal /= normal.two_norm();
      return normal;
    }

    //! Return the bound element stored as copy in the \ref bind function.
    const Element& localContext() const
    {
      return element_;
    }

    //! Not implemented.
    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:
    std::optional<Geometry> geometry_;
    Element element_;
  };

public:
  //! Construct the FaceNormalGridFunction.
  FaceNormalGridFunction(const GridView& gridView) :
    entitySet_(gridView)
  {}

  //! Not implemented.
  Range operator()(const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  //! Not implemented.
  friend typename Traits::DerivativeInterface derivative(const FaceNormalGridFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  //! Return a local-function associated to FaceNormalGridFunction.
  friend LocalFunction localFunction(const FaceNormalGridFunction& t)
  {
    return LocalFunction{};
  }

  //! Return the stored GridViewEntitySet.
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:
  EntitySet entitySet_;
};



} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_FACENORMALGRIDFUNCTION_HH
