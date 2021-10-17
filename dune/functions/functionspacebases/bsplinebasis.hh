// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH

/** \file
 * \brief The B-spline global function space basis
 */

#include <array>
#include <numeric>

/** \todo Don't use this matrix */
#include <dune/common/dynmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/common/diagonalmatrix.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/geometry/type.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

namespace Dune
{
namespace Functions {

// A maze of dependencies between the different parts of this.  We need a few forward declarations
template<typename GV, typename R>
class BSplineLocalFiniteElement;

template<typename GV>
class BSplinePreBasis;


/** \brief LocalBasis class in the sense of dune-localfunctions, presenting the restriction
 * of a B-spline patch to a knot span
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV Grid view that the basis is defined on
 * \tparam R Number type used for spline function values
 */
template<class GV, class R>
class BSplineLocalBasis
{
  friend class BSplineLocalFiniteElement<GV,R>;

  typedef typename GV::ctype D;
  enum {dim = GV::dimension};
public:

  //! \brief export type traits for function signature
  typedef LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,
  FieldMatrix<R,1,dim> > Traits;

  /** \brief Constructor with a given B-spline patch
   *
   * The patch object does all the work.
   */
  BSplineLocalBasis(const BSplinePreBasis<GV>& preBasis,
                    const BSplineLocalFiniteElement<GV,R>& lFE)
  : preBasis_(preBasis),
    lFE_(lFE)
  {}

  /** \brief Evaluate all shape functions
   * \param in Coordinates where to evaluate the functions, in local coordinates of the current knot span
   */
  void evaluateFunction (const FieldVector<D,dim>& in,
                         std::vector<FieldVector<R,1> >& out) const
  {
    FieldVector<D,dim> globalIn = offset_;
    scaling_.umv(in,globalIn);

    preBasis_.evaluateFunction(globalIn, out, lFE_.currentKnotSpan_);
  }

  /** \brief Evaluate Jacobian of all shape functions
   * \param in Coordinates where to evaluate the Jacobian, in local coordinates of the current knot span
   */
  void evaluateJacobian (const FieldVector<D,dim>& in,
                         std::vector<FieldMatrix<D,1,dim> >& out) const
  {
    FieldVector<D,dim> globalIn = offset_;
    scaling_.umv(in,globalIn);

    preBasis_.evaluateJacobian(globalIn, out, lFE_.currentKnotSpan_);

    for (size_t i=0; i<out.size(); i++)
      for (int j=0; j<dim; j++)
        out[i][0][j] *= scaling_[j][j];
  }

  //! \brief Evaluate all shape functions and derivatives of any order
  template<size_t k>
  inline void evaluate (const typename std::array<int,k>& directions,
                        const typename Traits::DomainType& in,
                        std::vector<typename Traits::RangeType>& out) const
  {
    switch(k)
    {
    case 0:
      evaluateFunction(in, out);
      break;
    case 1:
    {
      FieldVector<D,dim> globalIn = offset_;
      scaling_.umv(in,globalIn);

      preBasis_.evaluate(directions, globalIn, out, lFE_.currentKnotSpan_);

      for (size_t i=0; i<out.size(); i++)
        out[i][0] *= scaling_[directions[0]][directions[0]];
      break;
    }
    case 2:
    {
      FieldVector<D,dim> globalIn = offset_;
      scaling_.umv(in,globalIn);

      preBasis_.evaluate(directions, globalIn, out, lFE_.currentKnotSpan_);

      for (size_t i=0; i<out.size(); i++)
        out[i][0] *= scaling_[directions[0]][directions[0]]*scaling_[directions[1]][directions[1]];
      break;
    }
    default:
      DUNE_THROW(NotImplemented, "B-Spline derivatives of order " << k << " not implemented yet!");
    }
  }

  /** \brief Polynomial order of the shape functions
   *
   * Unfortunately, the general interface of the LocalBasis class mandates that the 'order' method
   * takes no arguments, and returns a single integer.  It therefore cannot reflect that fact that
   * a B-spline basis function can easily have different orders in the different coordinate directions.
   * We therefore take the conservative choice and return the maximum over the orders of all directions.
   */
  unsigned int order () const
  {
    return *std::max_element(preBasis_.order_.begin(), preBasis_.order_.end());
  }

  /** \brief Return the number of basis functions on the current knot span
   */
  std::size_t size() const
  {
    return lFE_.size();
  }

private:
  const BSplinePreBasis<GV>& preBasis_;

  const BSplineLocalFiniteElement<GV,R>& lFE_;

  // Coordinates in a single knot span differ from coordinates on the B-spline patch
  // by an affine transformation.  This transformation is stored in offset_ and scaling_.
  FieldVector<D,dim>    offset_;
  DiagonalMatrix<D,dim> scaling_;
};

/** \brief Attaches a shape function to an entity
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * The attachment uses the same order as for Qk elements.  This does *not* provide sufficient information
 * to compute global indices for the shape functions.  However, it does allow to find all degrees of freedom
 * that belong to the grid boundary, if the knot vector is open.
 *
 * \note Currently only implemented for 1d and 2d grids.  For higher dimensions you can still use
 *   the BSplineBasis, but you won't be able to find the degrees of freedom on the grid boundary.
 *
 * \tparam dim Dimension of the reference cube
 */
template<int dim>
class BSplineLocalCoefficients
{
  // Return i as a d-digit number in the (k+1)-nary system
  std::array<unsigned int,dim> multiindex (unsigned int i) const
  {
    std::array<unsigned int,dim> alpha;
    for (int j=0; j<dim; j++)
    {
      alpha[j] = i % sizes_[j];
      i = i/sizes_[j];
    }
    return alpha;
  }

  /** \brief Set the 'subentity' field for each dof for a 1d element */
  void setup1d(std::vector<unsigned int>& subEntity)
  {
    if (sizes_[0]==1)
    {
      subEntity[0] = 0;
      return;
    }

    /* edge and vertex numbering
       0----0----1
     */
    unsigned lastIndex=0;
    subEntity[lastIndex++] = 0;                 // corner 0
    for (unsigned i = 0; i < sizes_[0] - 2; ++i)
      subEntity[lastIndex++] = 0;               // inner dofs of element (0)

    subEntity[lastIndex++] = 1;                 // corner 1

    assert(size()==lastIndex);
  }

  void setup2d(std::vector<unsigned int>& subEntity)
  {
    unsigned lastIndex=0;

    // LocalKey: entity number , entity codim, dof indices within each entity
    /* edge and vertex numbering
       2----3----3
       |         |
       |         |
       0         1
       |         |
       |         |
       0----2----1
     */

    // lower edge (2)
    subEntity[lastIndex++] = 0;                 // corner 0
    for (unsigned i = 0; i < sizes_[0]-2; ++i)
      subEntity[lastIndex++] = 2;           // inner dofs of lower edge (2)

    subEntity[lastIndex++] = 1;                 // corner 1

    // iterate from bottom to top over inner edge dofs
    for (unsigned e = 0; e < sizes_[1]-2; ++e)
    {
      subEntity[lastIndex++] = 0;                   // left edge (0)
      for (unsigned i = 0; i < sizes_[0]-2; ++i)
        subEntity[lastIndex++] = 0;                     // face dofs
      subEntity[lastIndex++] = 1;                   // right edge (1)
    }

    // upper edge (3)
    subEntity[lastIndex++] = 2;                 // corner 2
    for (unsigned i = 0; i < sizes_[0]-2; ++i)
      subEntity[lastIndex++] = 3;                   // inner dofs of upper edge (3)

    subEntity[lastIndex++] = 3;                 // corner 3

    assert(size()==lastIndex);
  }


public:
  void init(const std::array<unsigned,dim>& sizes)
  {
    sizes_ = sizes;

    li_.resize(size());

    // Set up array of codimension-per-dof-number
    std::vector<unsigned int> codim(li_.size());

    for (std::size_t i=0; i<codim.size(); i++)
    {
      codim[i] = 0;
      // Codimension gets increased by 1 for each coordinate direction
      // where dof is on boundary
      std::array<unsigned int,dim> mIdx = multiindex(i);
      for (int j=0; j<dim; j++)
        if (mIdx[j]==0 or mIdx[j]==sizes[j]-1)
          codim[i]++;
    }

    // Set up index vector (the index of the dof in the set of dofs of a given subentity)
    // Algorithm: the 'index' has the same ordering as the dof number 'i'.
    // To make it consecutive we interpret 'i' in the (k+1)-adic system, omit all digits
    // that correspond to axes where the dof is on the element boundary, and transform the
    // rest to the (k-1)-adic system.
    std::vector<unsigned int> index(size());

    for (std::size_t i=0; i<index.size(); i++)
    {
      index[i] = 0;

      std::array<unsigned int,dim> mIdx = multiindex(i);

      for (int j=dim-1; j>=0; j--)
        if (mIdx[j]>0 and mIdx[j]<sizes[j]-1)
          index[i] = (sizes[j]-1)*index[i] + (mIdx[j]-1);
    }

    // Set up entity and dof numbers for each (supported) dimension separately
    std::vector<unsigned int> subEntity(li_.size());

    if (subEntity.size() > 0)
    {
      if (dim==1) {

        setup1d(subEntity);

      } else if (dim==2 and sizes_[0]>1 and sizes_[1]>1) {

        setup2d(subEntity);

      }
    }

    for (size_t i=0; i<li_.size(); i++)
      li_[i] = LocalKey(subEntity[i], codim[i], index[i]);
  }

  //! number of coefficients
  std::size_t size () const
  {
    return std::accumulate(sizes_.begin(), sizes_.end(), 1, std::multiplies<unsigned int>());
  }

  //! get i'th index
  const LocalKey& localKey (std::size_t i) const
  {
    return li_[i];
  }

private:

  // Number of shape functions on this element per coordinate direction
  std::array<unsigned, dim> sizes_;

  std::vector<LocalKey> li_;
};

/** \brief Local interpolation in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
template<int dim, class LB>
class BSplineLocalInterpolation
{
public:
  //! \brief Local interpolation of a function
  template<typename F, typename C>
  void interpolate (const F& f, std::vector<C>& out) const
  {
    DUNE_THROW(NotImplemented, "BSplineLocalInterpolation::interpolate");
  }

};

/** \brief LocalFiniteElement in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * This class ties together the implementation classes BSplineLocalBasis, BSplineLocalCoefficients, and BSplineLocalInterpolation
 *
 * \tparam D Number type used for domain coordinates
 * \tparam R Number type used for spline function values
 */
template<class GV, class R>
class BSplineLocalFiniteElement
{
  typedef typename GV::ctype D;
  enum {dim = GV::dimension};
  friend class BSplineLocalBasis<GV,R>;
public:

  /** \brief Export various types related to this LocalFiniteElement
   */
  typedef LocalFiniteElementTraits<BSplineLocalBasis<GV,R>,
  BSplineLocalCoefficients<dim>,
  BSplineLocalInterpolation<dim,BSplineLocalBasis<GV,R> > > Traits;

  /** \brief Constructor with a given B-spline basis
   */
  BSplineLocalFiniteElement(const BSplinePreBasis<GV>& preBasis)
  : preBasis_(preBasis),
    localBasis_(preBasis,*this)
  {}

  /** \brief Copy constructor
   */
  BSplineLocalFiniteElement(const BSplineLocalFiniteElement& other)
  : preBasis_(other.preBasis_),
    localBasis_(preBasis_,*this)
  {}

  /** \brief Bind LocalFiniteElement to a specific knot span of the spline patch
   *
   * Elements are the non-empty knot spans, here we do the renumbering
   *
   * \param ijk Integer coordinates in the tensor product patch
   */
  void bind(const std::array<unsigned,dim>& elementIdx)
  {
    /* \todo In the long run we need to precompute a table for this */
    for (size_t i=0; i<elementIdx.size(); i++)
    {
      currentKnotSpan_[i] = 0;

      // Skip over degenerate knot spans
      while (preBasis_.knotVectors_[i][currentKnotSpan_[i]+1] < preBasis_.knotVectors_[i][currentKnotSpan_[i]]+1e-8)
        currentKnotSpan_[i]++;

      for (size_t j=0; j<elementIdx[i]; j++)
      {
        currentKnotSpan_[i]++;

        // Skip over degenerate knot spans
        while (preBasis_.knotVectors_[i][currentKnotSpan_[i]+1] < preBasis_.knotVectors_[i][currentKnotSpan_[i]]+1e-8)
          currentKnotSpan_[i]++;
      }

      // Compute the geometric transformation from knotspan-local to global coordinates
      localBasis_.offset_[i] = preBasis_.knotVectors_[i][currentKnotSpan_[i]];
      localBasis_.scaling_[i][i] = preBasis_.knotVectors_[i][currentKnotSpan_[i]+1] - preBasis_.knotVectors_[i][currentKnotSpan_[i]];
    }

    // Set up the LocalCoefficients object
    std::array<unsigned int, dim> sizes;
    for (size_t i=0; i<dim; i++)
      sizes[i] = size(i);
    localCoefficients_.init(sizes);
  }

  /** \brief Hand out a LocalBasis object */
  const BSplineLocalBasis<GV,R>& localBasis() const
  {
    return localBasis_;
  }

  /** \brief Hand out a LocalCoefficients object */
  const BSplineLocalCoefficients<dim>& localCoefficients() const
  {
    return localCoefficients_;
  }

  /** \brief Hand out a LocalInterpolation object */
  const BSplineLocalInterpolation<dim,BSplineLocalBasis<GV,R> >& localInterpolation() const
  {
    return localInterpolation_;
  }

  /** \brief Number of shape functions in this finite element */
  unsigned size () const
  {
    std::size_t r = 1;
    for (int i=0; i<dim; i++)
      r *= size(i);
    return r;
  }

  /** \brief Return the reference element that the local finite element is defined on (here, a hypercube)
   */
  GeometryType type () const
  {
    return GeometryTypes::cube(dim);
  }

//private:

  /** \brief Number of degrees of freedom for one coordinate direction */
  unsigned int size(int i) const
  {
    const auto& order = preBasis_.order_;
    unsigned int r = order[i]+1;   // The 'normal' value
    if (currentKnotSpan_[i]<order[i])   // Less near the left end of the knot vector
      r -= (order[i] - currentKnotSpan_[i]);
    if ( order[i] > (preBasis_.knotVectors_[i].size() - currentKnotSpan_[i] - 2) )
      r -= order[i] - (preBasis_.knotVectors_[i].size() - currentKnotSpan_[i] - 2);
    return r;
  }

  const BSplinePreBasis<GV>& preBasis_;

  BSplineLocalBasis<GV,R> localBasis_;
  BSplineLocalCoefficients<dim> localCoefficients_;
  BSplineLocalInterpolation<dim,BSplineLocalBasis<GV,R> > localInterpolation_;

  // The knot span we are bound to
  std::array<unsigned,dim> currentKnotSpan_;
};


template<typename GV>
class BSplineNode;

/** \brief Pre-basis for B-spline basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 *
 * The BSplinePreBasis can be used to embed a BSplineBasis
 * in a larger basis for the construction of product spaces.
 */
template<typename GV>
class BSplinePreBasis
{
  static const int dim = GV::dimension;

  /** \brief Simple dim-dimensional multi-index class */
  class MultiDigitCounter
  {
  public:

    /** \brief Constructs a new multi-index, and sets all digits to zero
     *  \param limits Number of different digit values for each digit, i.e., digit i counts from 0 to limits[i]-1
     */
    MultiDigitCounter(const std::array<unsigned int,dim>& limits)
    : limits_(limits)
    {
      std::fill(counter_.begin(), counter_.end(), 0);
    }

    /** \brief Increment the multi-index */
    MultiDigitCounter& operator++()
    {
      for (int i=0; i<dim; i++)
      {
        ++counter_[i];

        // no overflow?
        if (counter_[i] < limits_[i])
          break;

        counter_[i] = 0;
      }
      return *this;
    }

    /** \brief Access the i-th digit of the multi-index */
    const unsigned int& operator[](int i) const
    {
      return counter_[i];
    }

    /** \brief How many times can you increment this multi-index before it overflows? */
    unsigned int cycle() const
    {
      unsigned int r = 1;
      for (int i=0; i<dim; i++)
        r *= limits_[i];
      return r;
    }

  private:

    /** \brief The number of different digit values for each place */
    const std::array<unsigned int,dim> limits_;

    /** \brief The values of the multi-index.  Each array entry is one digit */
    std::array<unsigned int,dim> counter_;

  };

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  using Node = BSplineNode<GV>;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  // Type used for function values
  using R = double;

  /** \brief Construct a B-spline basis for a given grid view and set of knot vectors
   *
   * The grid *must* match the knot vectors, i.e.:
   *  - The grid must be structured and Cartesian, and have cube elements only
   *  - The number of elements in each direction must match the number of knot spans in that direction
   *  - In fact, the element spacing in any direction must match the knot spacing in that direction
   *    (disregarding knot multiplicities)
   *  - When ordering the grid elements according to their indices, the resulting order must
   *    be lexicographical, with the x-index increasing fastest.
   *
   * Unfortunately, not all of these conditions can be checked for automatically.
   *
   * \param knotVector A single knot vector, which will be used for all coordinate directions
   * \param order B-spline order, will be used for all coordinate directions
   * \param makeOpen If this is true, then knots are prepended and appended to the knot vector to make the knot vector 'open'.
   *        i.e., start and end with 'order+1' identical knots.  Basis functions from such knot vectors are interpolatory at
   *        the end of the parameter interval.
   */
  BSplinePreBasis(const GridView& gridView,
                     const std::vector<double>& knotVector,
                     unsigned int order,
                     bool makeOpen = true)
  : gridView_(gridView)
  {
    // \todo Detection of duplicate knots
    std::fill(elements_.begin(), elements_.end(), knotVector.size()-1);

    // Mediocre sanity check: we don't know the number of grid elements in each direction.
    // but at least we know the total number of elements.
    assert( std::accumulate(elements_.begin(), elements_.end(), 1, std::multiplies<unsigned>()) == gridView_.size(0) );

    for (int i=0; i<dim; i++)
    {
      // Prepend the correct number of additional knots to open the knot vector
      //! \todo maybe test whether the knot vector is already open?
      if (makeOpen)
        for (unsigned int j=0; j<order; j++)
          knotVectors_[i].push_back(knotVector[0]);

      knotVectors_[i].insert(knotVectors_[i].end(), knotVector.begin(), knotVector.end());

      if (makeOpen)
        for (unsigned int j=0; j<order; j++)
          knotVectors_[i].push_back(knotVector.back());
    }

    std::fill(order_.begin(), order_.end(), order);
  }

  /** \brief Construct a B-spline basis for a given grid view with uniform knot vectors
   *
   * The grid *must* match the knot vectors, i.e.:
   *  - The grid must be structured and Cartesian, and have cube elements only
   *  - Bounding box and number of elements of the grid must match the corresponding arguments
   *    given to this constructor.
   *  - The element spacing must be uniform
   *  - When ordering the grid elements according to their indices, the resulting order must
   *    be lexicographical, with the x-index increasing fastest.
   *
   * Unfortunately, not all of these conditions can be checked for automatically.
   *
   * \param gridView The grid we are defining the basis on
   * \param lowerLeft Lower left corner of the structured grid
   * \param upperRight Upper right corner of the structured grid
   * \param elements Number of elements in each coordinate direction
   * \param order B-spline order, will be used for all coordinate directions
   * \param makeOpen If this is true, then knots are prepended and appended to the knot vector to make the knot vector 'open'.
   *        i.e., start and end with 'order+1' identical knots.  Basis functions from such knot vectors are interpolatory at
   *        the end of the parameter interval.
   */
  BSplinePreBasis(const GridView& gridView,
                     const FieldVector<double,dim>& lowerLeft,
                     const FieldVector<double,dim>& upperRight,
                     const std::array<unsigned int,dim>& elements,
                     unsigned int order,
                     bool makeOpen = true)
  : elements_(elements),
    gridView_(gridView)
  {
    // Mediocre sanity check: we don't know the number of grid elements in each direction.
    // but at least we know the total number of elements.
    assert( std::accumulate(elements_.begin(), elements_.end(), 1, std::multiplies<unsigned>()) == gridView_.size(0) );

    for (int i=0; i<dim; i++)
    {
      // Prepend the correct number of additional knots to open the knot vector
      //! \todo maybe test whether the knot vector is already open?
      if (makeOpen)
        for (unsigned int j=0; j<order; j++)
          knotVectors_[i].push_back(lowerLeft[i]);

      // Construct the actual knot vector
      for (size_t j=0; j<elements[i]+1; j++)
        knotVectors_[i].push_back(lowerLeft[i] + j*(upperRight[i]-lowerLeft[i]) / elements[i]);

      if (makeOpen)
        for (unsigned int j=0; j<order; j++)
          knotVectors_[i].push_back(upperRight[i]);
    }

    std::fill(order_.begin(), order_.end(), order);
  }

  //! Initialize the global indices
  void initializeIndices()
  {}

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
    gridView_ = gv;
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{this};
  }

  // Ideally this method should be implemented as
  //
  //   template<class SizePrefix>
  //   size_type size(const SizePrefix& prefix) const
  //
  // But leads to ambiguity with the other size method:
  //
  //   unsigned int size (size_t d) const
  //
  // Once the latter is removed, this implementation should be changed.

  //! Return number of possible values for next position in multi index
  template<class ST, int i>
  size_type size(const Dune::ReservedVector<ST, i>& prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    size_type result = 1;
    for (int i=0; i<dim; i++)
      result *= order_[i]+1;
    return result;
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it) const
  {
    // Local degrees of freedom are arranged in a lattice.
    // We need the lattice dimensions to be able to compute lattice coordinates from a local index
    std::array<unsigned int, dim> localSizes;
    for (int i=0; i<dim; i++)
      localSizes[i] = node.finiteElement().size(i);
    for (size_type i = 0, end = node.size() ; i < end ; ++i, ++it)
      {
        std::array<unsigned int,dim> localIJK = getIJK(i, localSizes);

        const auto currentKnotSpan = node.finiteElement().currentKnotSpan_;
        const auto order = order_;

        std::array<unsigned int,dim> globalIJK;
        for (int i=0; i<dim; i++)
          globalIJK[i] = std::max((int)currentKnotSpan[i] - (int)order[i], 0) + localIJK[i];  // needs to be a signed type!

        // Make one global flat index from the globalIJK tuple
        size_type globalIdx = globalIJK[dim-1];

        for (int i=dim-2; i>=0; i--)
          globalIdx = globalIdx * size(i) + globalIJK[i];

        *it = {{globalIdx}};
      }
    return it;
  }

  //! \brief Total number of B-spline basis functions
  unsigned int size () const
  {
    unsigned int result = 1;
    for (size_t i=0; i<dim; i++)
      result *= size(i);
    return result;
  }

  //! \brief Number of shape functions in one direction
  unsigned int size (size_t d) const
  {
    return knotVectors_[d].size() - order_[d] - 1;
  }

  /** \brief Evaluate all B-spline basis functions at a given point
   */
  void evaluateFunction (const FieldVector<typename GV::ctype,dim>& in,
                         std::vector<FieldVector<R,1> >& out,
                         const std::array<unsigned,dim>& currentKnotSpan) const
  {
    // Evaluate
    std::array<std::vector<R>, dim> oneDValues;

    for (size_t i=0; i<dim; i++)
      evaluateFunction(in[i], oneDValues[i], knotVectors_[i], order_[i], currentKnotSpan[i]);

    std::array<unsigned int, dim> limits;
    for (int i=0; i<dim; i++)
      limits[i] = oneDValues[i].size();

    MultiDigitCounter ijkCounter(limits);

    out.resize(ijkCounter.cycle());

    for (size_t i=0; i<out.size(); i++, ++ijkCounter)
    {
      out[i] = R(1.0);
      for (size_t j=0; j<dim; j++)
        out[i] *= oneDValues[j][ijkCounter[j]];
    }
  }

  /** \brief Evaluate Jacobian of all B-spline basis functions
   *
   * In theory this is easy: just look up the formula in a B-spline text of your choice.
   * The challenge is compute only the values needed for the current knot span.
   */
  void evaluateJacobian (const FieldVector<typename GV::ctype,dim>& in,
                         std::vector<FieldMatrix<R,1,dim> >& out,
                         const std::array<unsigned,dim>& currentKnotSpan) const
  {
    // How many shape functions to we have in each coordinate direction?
    std::array<unsigned int, dim> limits;
    for (int i=0; i<dim; i++)
    {
      limits[i] = order_[i]+1;  // The 'standard' value away from the boundaries of the knot vector
      if (currentKnotSpan[i]<order_[i])
        limits[i] -= (order_[i] - currentKnotSpan[i]);
      if ( order_[i] > (knotVectors_[i].size() - currentKnotSpan[i] - 2) )
        limits[i] -= order_[i] - (knotVectors_[i].size() - currentKnotSpan[i] - 2);
    }

    // The lowest knot spans that we need values from
    std::array<unsigned int, dim> offset;
    for (int i=0; i<dim; i++)
      offset[i] = std::max((int)(currentKnotSpan[i] - order_[i]),0);

    // Evaluate 1d function values (needed for the product rule)
    std::array<std::vector<R>, dim> oneDValues;

    // Evaluate 1d function values of one order lower (needed for the derivative formula)
    std::array<std::vector<R>, dim> lowOrderOneDValues;

    std::array<DynamicMatrix<R>, dim> values;

    for (size_t i=0; i<dim; i++)
    {
      evaluateFunctionFull(in[i], values[i], knotVectors_[i], order_[i], currentKnotSpan[i]);
      oneDValues[i].resize(knotVectors_[i].size()-order_[i]-1);
      for (size_t j=0; j<oneDValues[i].size(); j++)
        oneDValues[i][j] = values[i][order_[i]][j];

      if (order_[i]!=0)
      {
        lowOrderOneDValues[i].resize(knotVectors_[i].size()-(order_[i]-1)-1);
        for (size_t j=0; j<lowOrderOneDValues[i].size(); j++)
          lowOrderOneDValues[i][j] = values[i][order_[i]-1][j];
      }
    }


    // Evaluate 1d function derivatives
    std::array<std::vector<R>, dim> oneDDerivatives;
    for (size_t i=0; i<dim; i++)
    {
      oneDDerivatives[i].resize(limits[i]);

      if (order_[i]==0)  // order-zero functions are piecewise constant, hence all derivatives are zero
        std::fill(oneDDerivatives[i].begin(), oneDDerivatives[i].end(), R(0.0));
      else
      {
        for (size_t j=offset[i]; j<offset[i]+limits[i]; j++)
        {
          R derivativeAddend1 = lowOrderOneDValues[i][j] / (knotVectors_[i][j+order_[i]]-knotVectors_[i][j]);
          R derivativeAddend2 = lowOrderOneDValues[i][j+1] / (knotVectors_[i][j+order_[i]+1]-knotVectors_[i][j+1]);
          // The two previous terms may evaluate as 0/0.  This is to be interpreted as 0.
          if (std::isnan(derivativeAddend1))
            derivativeAddend1 = 0;
          if (std::isnan(derivativeAddend2))
            derivativeAddend2 = 0;
          oneDDerivatives[i][j-offset[i]] = order_[i] * ( derivativeAddend1 - derivativeAddend2 );
        }
      }
    }

    // Working towards computing only the parts that we really need:
    // Let's copy them out into a separate array
    std::array<std::vector<R>, dim> oneDValuesShort;

    for (int i=0; i<dim; i++)
    {
      oneDValuesShort[i].resize(limits[i]);

      for (size_t j=0; j<limits[i]; j++)
        oneDValuesShort[i][j]      = oneDValues[i][offset[i] + j];
    }



    // Set up a multi-index to go from consecutive indices to integer coordinates
    MultiDigitCounter ijkCounter(limits);

    out.resize(ijkCounter.cycle());

    // Complete Jacobian is given by the product rule
    for (size_t i=0; i<out.size(); i++, ++ijkCounter)
      for (int j=0; j<dim; j++)
      {
        out[i][0][j] = 1.0;
        for (int k=0; k<dim; k++)
          out[i][0][j] *= (j==k) ? oneDDerivatives[k][ijkCounter[k]]
                                 : oneDValuesShort[k][ijkCounter[k]];
      }

  }

  //! \brief Evaluate Derivatives of all B-spline basis functions
  template <size_type k>
  void evaluate(const typename std::array<int,k>& directions,
                const FieldVector<typename GV::ctype,dim>& in,
                std::vector<FieldVector<R,1> >& out,
                const std::array<unsigned,dim>& currentKnotSpan) const
  {
    if (k != 1 && k != 2)
      DUNE_THROW(RangeError, "Differentiation order greater than 2 is not supported!");

    // Evaluate 1d function values (needed for the product rule)
    std::array<std::vector<R>, dim> oneDValues;
    std::array<std::vector<R>, dim> oneDDerivatives;
    std::array<std::vector<R>, dim> oneDSecondDerivatives;

    // Evaluate 1d function derivatives
    if (k==1)
      for (size_t i=0; i<dim; i++)
        evaluateAll(in[i], oneDValues[i], true, oneDDerivatives[i], false, oneDSecondDerivatives[i], knotVectors_[i], order_[i], currentKnotSpan[i]);
    else
      for (size_t i=0; i<dim; i++)
        evaluateAll(in[i], oneDValues[i], true, oneDDerivatives[i], true, oneDSecondDerivatives[i], knotVectors_[i], order_[i], currentKnotSpan[i]);

    // The lowest knot spans that we need values from
    std::array<unsigned int, dim> offset;
    for (int i=0; i<dim; i++)
      offset[i] = std::max((int)(currentKnotSpan[i] - order_[i]),0);

    // Set up a multi-index to go from consecutive indices to integer coordinates
    std::array<unsigned int, dim> limits;
    for (int i=0; i<dim; i++)
    {
      // In a proper implementation, the following line would do
      //limits[i] = oneDValues[i].size();
      limits[i] = order_[i]+1;  // The 'standard' value away from the boundaries of the knot vector
      if (currentKnotSpan[i]<order_[i])
        limits[i] -= (order_[i] - currentKnotSpan[i]);
      if ( order_[i] > (knotVectors_[i].size() - currentKnotSpan[i] - 2) )
        limits[i] -= order_[i] - (knotVectors_[i].size() - currentKnotSpan[i] - 2);
    }

    // Working towards computing only the parts that we really need:
    // Let's copy them out into a separate array
    std::array<std::vector<R>, dim> oneDValuesShort;

    for (int i=0; i<dim; i++)
    {
      oneDValuesShort[i].resize(limits[i]);

      for (size_t j=0; j<limits[i]; j++)
        oneDValuesShort[i][j]  = oneDValues[i][offset[i] + j];
    }


    MultiDigitCounter ijkCounter(limits);

    out.resize(ijkCounter.cycle());

    if (k == 1)
    {
      // Complete Jacobian is given by the product rule
      for (size_t i=0; i<out.size(); i++, ++ijkCounter)
      {
        out[i][0] = 1.0;
        for (int l=0; l<dim; l++)
          out[i][0] *= (directions[0]==l) ? oneDDerivatives[l][ijkCounter[l]]
                                          : oneDValuesShort[l][ijkCounter[l]];
      }
    }

    if (k == 2)
    {
      // Complete derivation by deriving the tensor product
      for (size_t i=0; i<out.size(); i++, ++ijkCounter)
      {
        out[i][0] = 1.0;
        for (int j=0; j<dim; j++)
        {
          if (directions[0] != directions[1]) //derivation in two different variables
            if (directions[0] == j || directions[1] == j) //the spline has to be derived (once) in this direction
              out[i][0] *= oneDDerivatives[j][ijkCounter[j]];
            else //no derivation in this direction
              out[i][0] *= oneDValuesShort[j][ijkCounter[j]];
          else //spline is derived two times in the same direction
            if (directions[0] == j) //the spline is derived two times in this direction
              out[i][0] *= oneDSecondDerivatives[j][ijkCounter[j]];
            else //no derivation in this direction
              out[i][0] *= oneDValuesShort[j][ijkCounter[j]];
        }
      }
    }
  }


  /** \brief Compute integer element coordinates from the element index
   * \warning This method makes strong assumptions about the grid, namely that it is
   *   structured, and that indices are given in a x-fastest fashion.
   */
  static std::array<unsigned int,dim> getIJK(typename GridView::IndexSet::IndexType idx, std::array<unsigned int,dim> elements)
  {
    std::array<unsigned,dim> result;
    for (int i=0; i<dim; i++)
    {
      result[i] = idx%elements[i];
      idx /= elements[i];
    }
    return result;
  }

  /** \brief Evaluate all one-dimensional B-spline functions for a given coordinate direction
   *
   * This implementations was based on the explanations in the book of
   * Cottrell, Hughes, Bazilevs, "Isogeometric Analysis"
   *
   * \param in Scalar(!) coordinate where to evaluate the functions
   * \param [out] out Vector containing the values of all B-spline functions at 'in'
   */
  static void evaluateFunction (const typename GV::ctype& in, std::vector<R>& out,
                                const std::vector<R>& knotVector,
                                unsigned int order,
                                unsigned int currentKnotSpan)
  {
    std::size_t outSize = order+1;  // The 'standard' value away from the boundaries of the knot vector
    if (currentKnotSpan<order)   // Less near the left end of the knot vector
      outSize -= (order - currentKnotSpan);
    if ( order > (knotVector.size() - currentKnotSpan - 2) )
      outSize -= order - (knotVector.size() - currentKnotSpan - 2);
    out.resize(outSize);

    // It's not really a matrix that is needed here, a plain 2d array would do
    DynamicMatrix<R> N(order+1, knotVector.size()-1);

    // The text books on splines use the following geometric condition here to fill the array N
    // (see for example Cottrell, Hughes, Bazilevs, Formula (2.1).  However, this condition
    // only works if splines are never evaluated exactly on the knots.
    //
    // for (size_t i=0; i<knotVector.size()-1; i++)
    //   N[0][i] = (knotVector[i] <= in) and (in < knotVector[i+1]);
    for (size_t i=0; i<knotVector.size()-1; i++)
      N[0][i] = (i == currentKnotSpan);

    for (size_t r=1; r<=order; r++)
      for (size_t i=0; i<knotVector.size()-r-1; i++)
      {
        R factor1 = ((knotVector[i+r] - knotVector[i]) > 1e-10)
        ? (in - knotVector[i]) / (knotVector[i+r] - knotVector[i])
        : 0;
        R factor2 = ((knotVector[i+r+1] - knotVector[i+1]) > 1e-10)
        ? (knotVector[i+r+1] - in) / (knotVector[i+r+1] - knotVector[i+1])
        : 0;
        N[r][i] = factor1 * N[r-1][i] + factor2 * N[r-1][i+1];
      }

    /** \todo We only hand out function values for those basis functions whose support overlaps
     *  the current knot span.  However, in the preceding loop we still computed _all_ values_.
     * This won't scale.
     */
    for (size_t i=0; i<out.size(); i++) {
      out[i] = N[order][std::max((int)(currentKnotSpan - order),0) + i];
    }
  }

  /** \brief Evaluate all one-dimensional B-spline functions for a given coordinate direction
   *
   * This implementations was based on the explanations in the book of
   * Cottrell, Hughes, Bazilevs, "Isogeometric Analysis"
   *
   * \todo This method is a hack!  I computes the derivatives of ALL B-splines, even the ones that
   * are zero on the current knot span.  I need it as an intermediate step to get the derivatives
   * working.  It will/must be removed as soon as possible.
   *
   * \param in Scalar(!) coordinate where to evaluate the functions
   * \param [out] out Vector containing the values of all B-spline functions at 'in'
   */
  static void evaluateFunctionFull(const typename GV::ctype& in,
                                   DynamicMatrix<R>& out,
                                   const std::vector<R>& knotVector,
                                   unsigned int order,
                                   unsigned int currentKnotSpan)
  {
    // It's not really a matrix that is needed here, a plain 2d array would do
    DynamicMatrix<R>& N = out;

    N.resize(order+1, knotVector.size()-1);

    // The text books on splines use the following geometric condition here to fill the array N
    // (see for example Cottrell, Hughes, Bazilevs, Formula (2.1).  However, this condition
    // only works if splines are never evaluated exactly on the knots.
    //
    // for (size_t i=0; i<knotVector.size()-1; i++)
    //   N[0][i] = (knotVector[i] <= in) and (in < knotVector[i+1]);
    for (size_t i=0; i<knotVector.size()-1; i++)
      N[0][i] = (i == currentKnotSpan);

    for (size_t r=1; r<=order; r++)
      for (size_t i=0; i<knotVector.size()-r-1; i++)
      {
        R factor1 = ((knotVector[i+r] - knotVector[i]) > 1e-10)
        ? (in - knotVector[i]) / (knotVector[i+r] - knotVector[i])
        : 0;
        R factor2 = ((knotVector[i+r+1] - knotVector[i+1]) > 1e-10)
        ? (knotVector[i+r+1] - in) / (knotVector[i+r+1] - knotVector[i+1])
        : 0;
        N[r][i] = factor1 * N[r-1][i] + factor2 * N[r-1][i+1];
      }
  }


  /** \brief Evaluate the second derivatives of all one-dimensional B-spline functions for a given coordinate direction
   *
   * \param in Scalar(!) coordinate where to evaluate the functions
   * \param enableEvaluations switches calculation of desired derivatives on
   * \param [out] out Vector containing the values of all B-spline derivatives at 'in'
   * \param [out] outJac Vector containing the first derivatives of all B-spline derivatives at 'in' (only if calculation was switched on by enableEvaluations)
   * \param [out] outHess Vector containing the second derivatives of all B-spline derivatives at 'in' (only if calculation was switched on by enableEvaluations)
   */

  static void evaluateAll(const typename GV::ctype& in,
                                   std::vector<R>& out,
                                   bool evaluateJacobian, std::vector<R>& outJac,
                                   bool evaluateHessian, std::vector<R>& outHess,
                                   const std::vector<R>& knotVector,
                                   unsigned int order,
                                   unsigned int currentKnotSpan)
  {
    // How many shape functions to we have in each coordinate direction?
    unsigned int limit;
    limit = order+1;  // The 'standard' value away from the boundaries of the knot vector
    if (currentKnotSpan<order)
      limit -= (order - currentKnotSpan);
    if ( order > (knotVector.size() - currentKnotSpan - 2) )
      limit -= order - (knotVector.size() - currentKnotSpan - 2);

    // The lowest knot spans that we need values from
    unsigned int offset;
    offset = std::max((int)(currentKnotSpan - order),0);

    // Evaluate 1d function values (needed for the product rule)
    DynamicMatrix<R> values;

    evaluateFunctionFull(in, values, knotVector, order, currentKnotSpan);

    out.resize(knotVector.size()-order-1);
    for (size_t j=0; j<out.size(); j++)
        out[j] = values[order][j];

    // Evaluate 1d function values of one order lower (needed for the derivative formula)
    std::vector<R> lowOrderOneDValues;

    if (order!=0)
    {
      lowOrderOneDValues.resize(knotVector.size()-(order-1)-1);
      for (size_t j=0; j<lowOrderOneDValues.size(); j++)
        lowOrderOneDValues[j] = values[order-1][j];
    }

    // Evaluate 1d function values of two order lower (needed for the (second) derivative formula)
    std::vector<R> lowOrderTwoDValues;

    if (order>1 && evaluateHessian)
    {
      lowOrderTwoDValues.resize(knotVector.size()-(order-2)-1);
      for (size_t j=0; j<lowOrderTwoDValues.size(); j++)
        lowOrderTwoDValues[j] = values[order-2][j];
    }

    // Evaluate 1d function derivatives
    if (evaluateJacobian)
    {
      outJac.resize(limit);

      if (order==0)  // order-zero functions are piecewise constant, hence all derivatives are zero
        std::fill(outJac.begin(), outJac.end(), R(0.0));
      else
      {
        for (size_t j=offset; j<offset+limit; j++)
        {
          R derivativeAddend1 = lowOrderOneDValues[j] / (knotVector[j+order]-knotVector[j]);
          R derivativeAddend2 = lowOrderOneDValues[j+1] / (knotVector[j+order+1]-knotVector[j+1]);
          // The two previous terms may evaluate as 0/0.  This is to be interpreted as 0.
          if (std::isnan(derivativeAddend1))
            derivativeAddend1 = 0;
          if (std::isnan(derivativeAddend2))
            derivativeAddend2 = 0;
          outJac[j-offset] = order * ( derivativeAddend1 - derivativeAddend2 );
        }
      }
    }

    // Evaluate 1d function second derivatives
    if (evaluateHessian)
    {
      outHess.resize(limit);

      if (order<2)  // order-zero functions are piecewise constant, hence all derivatives are zero
        std::fill(outHess.begin(), outHess.end(), R(0.0));
      else
      {
        for (size_t j=offset; j<offset+limit; j++)
        {
          assert(j+2 < lowOrderTwoDValues.size());
          R derivativeAddend1 = lowOrderTwoDValues[j] / (knotVector[j+order]-knotVector[j]) / (knotVector[j+order-1]-knotVector[j]);
          R derivativeAddend2 = lowOrderTwoDValues[j+1] / (knotVector[j+order]-knotVector[j]) / (knotVector[j+order]-knotVector[j+1]);
          R derivativeAddend3 = lowOrderTwoDValues[j+1] / (knotVector[j+order+1]-knotVector[j+1]) / (knotVector[j+order]-knotVector[j+1]);
          R derivativeAddend4 = lowOrderTwoDValues[j+2] / (knotVector[j+order+1]-knotVector[j+1]) / (knotVector[j+1+order]-knotVector[j+2]);
          // The two previous terms may evaluate as 0/0.  This is to be interpreted as 0.

          if (std::isnan(derivativeAddend1))
            derivativeAddend1 = 0;
          if (std::isnan(derivativeAddend2))
            derivativeAddend2 = 0;
          if (std::isnan(derivativeAddend3))
            derivativeAddend3 = 0;
          if (std::isnan(derivativeAddend4))
            derivativeAddend4 = 0;
          outHess[j-offset] = order * (order-1) * ( derivativeAddend1 - derivativeAddend2 -derivativeAddend3 + derivativeAddend4 );
        }
      }
    }
  }


  /** \brief Order of the B-spline for each space dimension */
  std::array<unsigned int, dim> order_;

  /** \brief The knot vectors, one for each space dimension */
  std::array<std::vector<double>, dim> knotVectors_;

  /** \brief Number of grid elements in the different coordinate directions */
  std::array<unsigned,dim> elements_;

  GridView gridView_;
};



template<typename GV>
class BSplineNode :
  public LeafBasisNode
{
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = BSplineLocalFiniteElement<GV,double>;

  BSplineNode(const BSplinePreBasis<GV>* preBasis) :
    preBasis_(preBasis),
    finiteElement_(*preBasis)
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return element_;
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
    element_ = e;
    auto elementIndex = preBasis_->gridView().indexSet().index(e);
    finiteElement_.bind(preBasis_->getIJK(elementIndex,preBasis_->elements_));
    this->setSize(finiteElement_.size());
  }

protected:

  const BSplinePreBasis<GV>* preBasis_;

  FiniteElement finiteElement_;
  Element element_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a B-spline pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 */
inline auto bSpline(const std::vector<double>& knotVector,
                    unsigned int order,
                    bool makeOpen = true)
{
  return [&knotVector, order, makeOpen](const auto& gridView) {
    return BSplinePreBasis<std::decay_t<decltype(gridView)>>(gridView, knotVector, order, makeOpen);
  };
}

} // end namespace BasisFactory

// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief A global B-spline basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using BSplineBasis = DefaultGlobalBasis<BSplinePreBasis<GV> >;


}   // namespace Functions

}   // namespace Dune

#endif   // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH
