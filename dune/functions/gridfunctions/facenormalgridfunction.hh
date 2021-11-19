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
 * \brief Grid function implementing the piecewise element normal
 *
 * This function computes implements the grid interface.
 * When evaluated at a point x inside of an element, it computes
 * the unit outward normal vector of the face closest to this
 * point.
 *
 * \ingroup FunctionImplementations
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

    void bind(const Element& element)
    {
      element_ = element;
      geometry_.emplace(element_.geometry());
    }

    void unbind()
    {}

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

    const Element& localContext() const
    {
      return element_;
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:
    std::optional<Geometry> geometry_;
    Element element_;
  };

public:

  FaceNormalGridFunction(const GridView& gridView) :
    entitySet_(gridView)
  {}

  Range operator()(const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend typename Traits::DerivativeInterface derivative(const FaceNormalGridFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend LocalFunction localFunction(const FaceNormalGridFunction& t)
  {
    return LocalFunction{};
  }

  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:
  EntitySet entitySet_;
};



} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_FACENORMALGRIDFUNCTION_HH
