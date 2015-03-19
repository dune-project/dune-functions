#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH

/** \file
 * \brief The B-spline global function space basis
 */

/** \todo Don't use this matrix */
#include <dune/common/dynmatrix.hh>

#include <dune/common/version.hh>

namespace Dune
{
namespace Functions {

// A maze of dependencies between the different parts of this.  We need lots of forward declarations
template<typename D, typename R, int dim>
class BSplineLocalFiniteElement;

template<typename GV>
class BSplineBasisLocalView;

template<typename GV>
class BSplineBasisLeafNode;

template<typename GV>
class BSplineIndexSet;

template <class D, class R, int dim>
class BSplineLocalBasis;

template<typename GV>
class BSplineBasis;

template<class D, class R, int dim>
class BSplinePatch;

/** \brief Maps local shape functions to global indices */
template<typename GV>
class BSplineLocalIndexSet
{
public:
  /** \brief Type used for sizes and local indices */
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef BSplineBasisLocalView<GV> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  BSplineLocalIndexSet(const BSplineIndexSet<GV> & indexSet)
  : indexSet_(indexSet)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const BSplineBasisLocalView<GV>& localView)
  {
    localView_ = &localView;
  }

  /** \brief Unbind the index set
   */
  void unbind()
  {
    localView_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return localView_->size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis (pair of multi-indices)
  MultiIndex index(size_type i) const
  {
    const auto currentKnotSpan = localView_->tree().finiteElement().currentKnotSpan_;
    const auto order = localView_->globalBasis_->patch_.order_;

    int offset = currentKnotSpan[0] - order[0];  // needs to be a signed type!
    offset = std::max(offset, 0);
    return { offset + i};
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

private:
  const BSplineBasisLocalView<GV>* localView_;

  const BSplineIndexSet<GV> indexSet_;
};

/** \brief Provides the size of the global basis, and hands out the local index sets */
template<typename GV>
class BSplineIndexSet
{
public:

  typedef BSplineLocalIndexSet<GV> LocalIndexSet;

  BSplineIndexSet(const BSplinePatch<typename GV::ctype, double, GV::dimension>* patch)
  : patch_(patch)
  {}

  /** \brief Total number of basis vectors in the basis */
  std::size_t size() const
  {
    return patch_->size();
  }

  /** \brief Provide a local index set, which hands out global indices for all shape functions of an element */
  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(*this);
  }

private:
  const BSplinePatch<typename GV::ctype, double, GV::dimension>* patch_;
};



/** \brief The restriction of a finite element basis to a single element */
template<typename GV>
class BSplineBasisLocalView
{
  // Grid dimension
  enum {dim = GV::dimension};

  // Needs the grid element
  friend class BSplineLocalIndexSet<GV>;

public:
  /** \brief The global FE basis that this is a view on */
  typedef BSplineBasis<GV> GlobalBasis;

  /** \brief The grid view of the global basis */
  typedef typename GlobalBasis::GridView GridView;

  /** \brief The type used for sizes */
  typedef typename GlobalBasis::size_type size_type;

  /** \brief Type used to number the global basis vectors
   *
   * In the case of mixed finite elements this really can be a multi-index, but for a standard
   * B-spline space this is only a single-digit multi-index, i.e., it is an integer.
   */
  typedef typename GlobalBasis::MultiIndex MultiIndex;

  /** \brief Type of the grid element we are bound to */
  typedef typename GridView::template Codim<0>::Entity Element;

  /** \brief Tree of local finite elements / local shape function sets
   *
   * In the case of a P2 space this tree consists of a single leaf only,
   * i.e., Tree is basically the type of the LocalFiniteElement
   */
  typedef BSplineBasisLeafNode<GV> Tree;

  /** \brief Construct local view for a given global finite element basis */
  BSplineBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(globalBasis)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = &e;
    tree_.bind(e);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   * And indeed, in the BSplineBasisView implementation this method does nothing.
   */
  void unbind()
  {}

  /** \brief Return the local ansatz tree associated to the bound entity
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /** \brief Number of degrees of freedom on this element
   */
  size_type size() const
  {
    return tree_.size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   */
  size_type maxSize() const
  {
    /** \todo This is not correct for non-open knot vectors */
    return size();
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  const GlobalBasis* globalBasis_;
  const Element* element_;
  Tree tree_;
};


/** \brief The finite element tree of this basis consists of a single node, and this is it */
template<typename GV>
class BSplineBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    BSplineLocalFiniteElement<typename GV::ctype,double,GV::dimension>,
    typename BSplineBasis<GV>::size_type>
{
  typedef BSplineBasis<GV> GlobalBasis;
  static const int dim = GV::dimension;

  typedef typename GV::template Codim<0>::Entity E;
  typedef BSplineLocalFiniteElement<typename GV::ctype,double,GV::dimension> FE;
  typedef typename GlobalBasis::size_type ST;
  typedef typename GlobalBasis::MultiIndex MI;

  typedef typename GlobalBasis::LocalView LocalView;

  friend LocalView;
  friend class BSplineLocalIndexSet<GV>;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;

  /** \brief Construct a leaf node for a given global B-spline basis */
  BSplineBasisLeafNode(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    finiteElement_(globalBasis->patch_),
    element_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const DUNE_FINAL
  {
    return finiteElement_;
  }

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type size() const DUNE_FINAL
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_.size();
  }

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const DUNE_FINAL
  {
    return i;
  }

private:
  /** \brief Bind to an element
   *
   * This involves in particular computing integer indices in a structured grid from the single
   * element index that dune-grid gives us.  Hence we have to make a few assumptions about the
   * grid itself, which is dangerous but cannot be helped.
   */
  void bind(const Element& e)
  {
    element_ = &e;

    auto elementIndex = globalBasis_->gridView().indexSet().index(e);
    finiteElement_.bind(globalBasis_->getIJK(elementIndex));
  }

  const GlobalBasis* globalBasis_;
  FiniteElement finiteElement_;
  const Element* element_;
};


/** \brief Tensor product of 1d B-Spline functions with arbitrary knot vectors
 *
 * \internal This is an internal implementation class, used to implement the BSplineBasis class.
 * \todo Try to merge this into the BSplineBasis class
 *
 * \tparam D Number type used for domain coordinates
 * \tparam R Number type used for spline function values
 * \tparam dim Dimension of the patch
 */
template<class D, class R, int dim>
class BSplinePatch
{
  template <class GV>
  friend class BSplineLocalIndexSet;
  friend class BSplineLocalFiniteElement<D,R,dim>;
  friend class BSplineLocalBasis<D,R,dim>;

  /** \brief Simple dim-dimensional multi-index class */
  class MultiIndex
  {
  public:

    /** \brief Constructs a new multi-index, and sets all digits to zero
     *  \param limits Number of different digit values for each digit, i.e., digit i counts from 0 to limits[i]-1
     */
    MultiIndex(const std::array<unsigned int,dim>& limits)
    : limits_(limits)
    {
      std::fill(counter_.begin(), counter_.end(), 0);
    }

    /** \brief Increment the multi-index */
    MultiIndex& operator++()
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

  /** \brief Constructor: same knot vector and order in all space directions
   * \param makeOpen If this is true, then knots are prepended and appended to the knot vector to make the knot vector 'open'.
   *        i.e., start and end with 'order+1' identical knots.  Basis functions from such knot vectors are interpolatory at
   *        the end of the parameter interval.
   */
  BSplinePatch(const std::vector<R>& knotVector,
               unsigned int order,
               bool makeOpen)
  {
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

  //! \brief Total number of B-spline basis functions
  unsigned int size () const
  {
    unsigned int result = 1;
    for (size_t i=0; i<dim; i++)
      result *= size(i);
    return result;
  }

  /** \brief Evaluate all B-spline basis functions at a given point
   */
  void evaluateFunction (const FieldVector<D,dim>& in,
                         std::vector<FieldVector<R,1> >& out,
                         const std::array<uint,dim>& currentKnotSpan) const
  {
    // Evaluate
    Dune::array<std::vector<R>, dim> oneDValues;

    for (size_t i=0; i<dim; i++)
      evaluateFunction(in[i], oneDValues[i], knotVectors_[i], order_[i], currentKnotSpan[i]);

    std::array<unsigned int, dim> limits;
    for (int i=0; i<dim; i++)
      limits[i] = oneDValues[i].size();

    MultiIndex ijkCounter(limits);

    out.resize(ijkCounter.cycle());

    for (size_t i=0; i<out.size(); i++, ++ijkCounter)
    {
      out[i] = R(1.0);
      for (size_t j=0; j<dim; j++)
        out[i] *= oneDValues[j][ijkCounter[j]];
    }
  }

  //! \brief Evaluate Jacobian of all B-spline basis functions
  void evaluateJacobian (const FieldVector<D,dim>& in,
                         std::vector<FieldMatrix<D,1,dim> >& out,
                         const std::array<uint,dim>& currentKnotSpan) const
  {
    // Evaluate 1d function values (needed for the product rule)
    Dune::array<std::vector<R>, dim> oneDValues;

    for (size_t i=0; i<dim; i++)
      evaluateFunctionFull(in[i], oneDValues[i], knotVectors_[i], order_[i], currentKnotSpan[i]);

    // Evaluate 1d function values of one order lower (needed for the derivative formula)
    Dune::array<std::vector<R>, dim> lowOrderOneDValues;

    /** \todo Calling evaluateFunction again here is a waste: the lower-order values have
     * already been computed during the first call to evaluateFunction.  Still, for the time
     * being I leave it as is to have more readable code. */
    for (size_t i=0; i<dim; i++)
      if (order_[i]!=0)
        evaluateFunctionFull(in[i], lowOrderOneDValues[i], knotVectors_[i], order_[i]-1, currentKnotSpan[i]);

    // Evaluate 1d function derivatives
    Dune::array<std::vector<R>, dim> oneDDerivatives;
    for (size_t i=0; i<dim; i++)
    {
      oneDDerivatives[i].resize(oneDValues[i].size());
      for (size_t j=0; j<oneDDerivatives[i].size(); j++)
      {
        if (order_[i]==0)  // order-zero functions are piecewise constant, hence all derivatives are zero
          std::fill(oneDDerivatives[i].begin(), oneDDerivatives[i].end(), R(0.0));
        else
          oneDDerivatives[i][j] = order_[i] * (lowOrderOneDValues[i][j] / (knotVectors_[i][j+order_[i]]-knotVectors_[i][j])
                                - lowOrderOneDValues[i][j+1] / (knotVectors_[i][j+order_[i]+1]-knotVectors_[i][j+1]) );

      }
    }

    // Set up a multi-index to go from consecutive indices to integer coordinates
    std::array<unsigned int, dim> limits;
    for (int i=0; i<dim; i++)
    {
      // In a proper implementation, the following line would do
      //limits[i] = oneDValues[i].size();
      limits[i] = order_[i]+1;  // The 'standard' value away from the boundaries of the knot vector
      limits[i] = std::min(limits[i], (unsigned)currentKnotSpan[i]+1);  // Less near the left end of the knot vector
      limits[i] = std::min(limits[i], (unsigned)knotVectors_[i].size() - currentKnotSpan[i]-1);  // Less near the right end of the knot vector
    }

    MultiIndex ijkCounter(limits);

    out.resize(ijkCounter.cycle());

    // Complete Jacobian is given by the product rule
    for (size_t i=0; i<out.size(); i++, ++ijkCounter)
      for (int j=0; j<dim; j++)
      {
        out[i][0][j] = 1.0;
        for (int k=0; k<dim; k++)
          out[i][0][j] *= (j==k) ? oneDDerivatives[k][std::max((int)(currentKnotSpan[k] - order_[k]),0) + ijkCounter[k]]
                                 : oneDValues[k][std::max((int)(currentKnotSpan[k] - order_[k]),0) + ijkCounter[k]];
      }

  }

  /** \brief Polynomial order of the shape functions
   * \todo Only implemented for 1d patches
   */
  unsigned int order () const
  {
    return order_;
  }

private:

  //! \brief Number of shape functions in one direction
  unsigned int size (size_t d) const
  {
    return knotVectors_[d].size() - order_[d] - 1;
  }


  /** \brief Evaluate all one-dimensional B-spline functions for a given coordinate direction
   *
   * This implementations was based on the explanations in the book of
   * Cottrell, Hughes, Bazilevs, "Isogeometric Analysis"
   *
   * \param in Scalar(!) coordinate where to evaluate the functions
   * \param [out] out Vector containing the values of all B-spline functions at 'in'
   */
  static void evaluateFunction (const D& in, std::vector<R>& out,
                                const std::vector<R>& knotVector,
                                unsigned int order,
                                unsigned int currentKnotSpan)
  {
    std::size_t outSize = order+1;  // The 'standard' value away from the boundaries of the knot vector
    outSize = std::min(outSize, (std::size_t)currentKnotSpan+1);  // Less near the left end of the knot vector
    outSize = std::min(outSize, knotVector.size() - currentKnotSpan-1);  // Less near the right end of the knot vector
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
  static void evaluateFunctionFull(const D& in, std::vector<R>& out,
                                   const std::vector<R>& knotVector,
                                   unsigned int order,
                                   unsigned int currentKnotSpan)
  {
    out.resize(knotVector.size()-order-1);

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

      for (size_t i=0; i<out.size(); i++) {
        out[i] = N[order][i];
      }
  }

  /** \brief Order of the B-spline for each space dimension */
  array<unsigned int, dim> order_;

  /** \brief The knot vectors, one for each space dimension */
  array<std::vector<R>, dim> knotVectors_;
};


/** \brief LocalBasis class in the sense of dune-localfunctions, presenting the restriction
 * of a B-spline patch to a knot span
 *
 * \tparam D Number type used for domain coordinates
 * \tparam R Number type used for spline function values
 * \tparam dim Dimension of the patch
 */
template<class D, class R, int dim>
class BSplineLocalBasis
{
  friend class BSplineLocalFiniteElement<D,R,dim>;
public:

  //! \brief export type traits for function signature
  typedef LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>,
  Dune::FieldMatrix<R,1,dim>, 2> Traits;

  /** \brief Constructor with a given B-spline patch
   *
   * The patch object does all the work.
   */
  BSplineLocalBasis(const BSplinePatch<D,R,dim>& patch,
                    const BSplineLocalFiniteElement<D,R,dim>& lFE)
  : patch_(patch),
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

    patch_.evaluateFunction(globalIn, out, lFE_.currentKnotSpan_);
  }

  /** \brief Evaluate Jacobian of all shape functions
   * \param in Coordinates where to evaluate the Jacobian, in local coordinates of the current knot span
   */
  void evaluateJacobian (const FieldVector<D,dim>& in,
                         std::vector<FieldMatrix<D,1,dim> >& out) const
  {
    FieldVector<D,dim> globalIn = offset_;
    scaling_.umv(in,globalIn);

    patch_.evaluateJacobian(globalIn, out, lFE_.currentKnotSpan_);

    for (size_t i=0; i<out.size(); i++)
      for (int j=0; j<dim; j++)
        out[i][0][j] *= scaling_[j][j];
  }

  //! \brief Evaluate all shape functions and derivatives of any order
  template<unsigned int k>
  inline void evaluate (const typename Dune::array<int,k>& directions,
                        const typename Traits::DomainType& in,
                        std::vector<typename Traits::RangeType>& out) const
  {
    if (k==0)
      evaluateFunction(in, out);
    else
      DUNE_THROW(NotImplemented, "B-Spline derivatives not implemented yet!");
  }

  /** \brief Polynomial order of the shape functions
   */
  unsigned int order () const
  {
    DUNE_THROW(NotImplemented, "!");
  }

  /** \brief Return the number of basis functions on the current knot span
   */
  std::size_t size() const
  {
    return lFE_.size();
  }

private:
  const BSplinePatch<D,R,dim>& patch_;

  const BSplineLocalFiniteElement<D,R,dim>& lFE_;

  // Coordinates in a single knot span differ from coordinates on the B-spline patch
  // by an affine transformation.  This transformation is stored in offset_ and scaling_.
  FieldVector<D,dim>    offset_;
  DiagonalMatrix<D,dim> scaling_;
};

/** \brief Local coefficients in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
 *
 *      \implements Dune::LocalCoefficientsVirtualImp
 */
template <class D, class R, int dim>
class BSplineLocalCoefficients
{
public:
  /** \brief Standard constructor
   * \todo Not implemented yet
   */
  BSplineLocalCoefficients (const BSplineLocalFiniteElement<D,R,dim>& lFE)
  : lFE_(lFE)
  {
    std::cout << "WARNING: LocalCoefficients array should be initialized here!" << std::endl;
  }

  /** \brief Number of coefficients
   * \todo Currently, the number of all basis functions on the entire patch is returned.
   *   This will include many that are simply constant zero on the current knot span.
   */
  std::size_t size () const
  {
    return lFE_.size();
  }

  //! get i'th index
  const LocalKey& localKey (std::size_t i) const
  {
    return li_[i];
  }

private:
  const BSplineLocalFiniteElement<D,R,dim>& lFE_;
  std::vector<LocalKey> li_;
};

/** \brief Local interpolation in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
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
 * This class ties together the implementation classes BSplineLocalBasis, BSplineLocalCoefficients, and BSplineLocalInterpolation
 *
 * \tparam D Number type used for domain coordinates
 * \tparam R Number type used for spline function values
 * \tparam dim Dimension of the patch
 */
template<class D, class R, int dim>
class BSplineLocalFiniteElement
{
  template <class GV>
  friend class BSplineLocalIndexSet;
  friend class BSplineLocalBasis<D,R,dim>;
public:

  /** \brief Export various types related to this LocalFiniteElement
   */
  typedef LocalFiniteElementTraits<BSplineLocalBasis<D,R,dim>,
  BSplineLocalCoefficients<D,R,dim>,
  BSplineLocalInterpolation<dim,BSplineLocalBasis<D,R,dim> > > Traits;

  /** \brief Constructor with a given B-spline patch
   */
  BSplineLocalFiniteElement(const BSplinePatch<D,R,dim>& patch)
  : patch_(patch),
    localBasis_(patch,*this),
    localCoefficients_(*this)
  {}

  /** \brief Bind LocalFiniteElement to a specific knot span of the spline patch
   *
   * Elements are the non-empty knot spans, here we do the renumbering
   *
   * \param ijk Integer coordinates in the tensor product patch
   */
  void bind(const std::array<uint,dim>& elementIdx)
  {
    /* \todo In the long run we need to precompute a table for this */
    for (size_t i=0; i<elementIdx.size(); i++)
    {
      currentKnotSpan_[i] = 0;

      // Skip over degenerate knot spans
      while (patch_.knotVectors_[i][currentKnotSpan_[i]+1] < patch_.knotVectors_[i][currentKnotSpan_[i]]+1e-8)
        currentKnotSpan_[i]++;

      for (size_t j=0; j<elementIdx[i]; j++)
      {
        currentKnotSpan_[i]++;

        // Skip over degenerate knot spans
        while (patch_.knotVectors_[i][currentKnotSpan_[i]+1] < patch_.knotVectors_[i][currentKnotSpan_[i]]+1e-8)
          currentKnotSpan_[i]++;
      }

      // Compute the geometric transformation from knotspan-local to global coordinates
      localBasis_.offset_[i] = patch_.knotVectors_[i][currentKnotSpan_[i]];
      localBasis_.scaling_[i][i] = patch_.knotVectors_[i][currentKnotSpan_[i]+1] - patch_.knotVectors_[i][currentKnotSpan_[i]];
    }
  }

  /** \brief Hand out a LocalBasis object */
  const BSplineLocalBasis<D,R,dim>& localBasis() const
  {
    return localBasis_;
  }

  /** \brief Hand out a LocalCoefficients object */
  const BSplineLocalCoefficients<D,R,dim>& localCoefficients() const
  {
    return localCoefficients_;
  }

  /** \brief Hand out a LocalInterpolation object */
  const BSplineLocalInterpolation<dim,BSplineLocalBasis<D,R,dim> >& localInterpolation() const
  {
    return localInterpolation_;
  }

  /** \brief Number of shape functions in this finite element */
  uint size () const
  {
    std::size_t r = 1;
    for (int i=0; i<dim; i++)
    {
      std::size_t oneDSize = patch_.order_[i]+1;   // The 'normal' value
      oneDSize = std::min(oneDSize, (std::size_t)currentKnotSpan_[i]+1);  // Less near the left end of the knot vector
      oneDSize = std::min(oneDSize, patch_.knotVectors_[i].size() - currentKnotSpan_[i]-1);  // Less near the right end of the knot vector
      r *= oneDSize;
    }
    return r;
  }

  /** \brief Return the reference element that the local finite element is defined on (here, a hypercube)
   */
  GeometryType type () const
  {
    return GeometryType(GeometryType::cube,dim);
  }

private:

  const BSplinePatch<D,R,dim>& patch_;
  BSplineLocalBasis<D,R,dim> localBasis_;
  BSplineLocalCoefficients<D,R,dim> localCoefficients_;
  BSplineLocalInterpolation<dim,BSplineLocalBasis<D,R,dim> > localInterpolation_;

  // The knot span we are bound to
  std::array<uint,dim> currentKnotSpan_;
};

/** \brief A B-spline function space basis on a tensor-product grid
 *
 * \tparam GV GridView, this must match the knot vectors describing the B-spline basis
 * \tparam RT Number type used for function values
 *
 * \todo Various features are not implemented yet:
 *  - No multiple knots in a knot vector
 *  - No sparsity; currently the implementation pretends that the support of any B-spline basis
 *    function covers the entire patch.
 */
template <class GV>
class BSplineBasis
: public GridViewFunctionSpaceBasis<GV,
                                    BSplineBasisLocalView<GV>,
                                    BSplineIndexSet<GV>,
                                    std::array<std::size_t, 1> >
{

  enum {dim = GV::dimension};

  friend class BSplineBasisLeafNode<GV>;
  friend class BSplineLocalIndexSet<GV>;

public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;

  /** \todo Do we really have to export this here? */
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef BSplineBasisLocalView<GV> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

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
  BSplineBasis(const GridView& gridView,
               const std::vector<double>& knotVector,
               unsigned int order,
               bool makeOpen = true)
  : patch_(knotVector, order, makeOpen),
    gridView_(gridView),
    indexSet_(&patch_)
  {
    // \todo Detection of duplicate knots
    std::fill(elements_.begin(), elements_.end(), knotVector.size()-1);

    // Mediocre sanity check: we don't know the number of grid elements in each direction.
    // but at least we now the total number of elements.
    assert( std::accumulate(elements_.begin(), elements_.end(), 1, std::multiplies<uint>()) == gridView_.size(0) );
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  BSplineIndexSet<GV> indexSet() const
  {
    return indexSet_;
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

protected:

  /** \brief Compute integer element coordinates from the element index
   * \warning This method makes strong assumptions about the grid, namely that it is
   *   structured, and that indices are given in a x-fastest fashion.
   */
  std::array<unsigned int,dim> getIJK(typename GridView::IndexSet::IndexType idx) const
  {
    std::array<uint,dim> result;
    for (int i=0; i<dim; i++)
    {
      result[i] = idx%elements_[i];
      idx /= elements_[i];
    }
    return result;
  }

  /** \brief The underlying B-spline basis patch */
  BSplinePatch<double,double,dim> patch_;

  /** \brief Number of grid elements in the different coordinate directions */
  std::array<uint,dim> elements_;

  const GridView gridView_;

  BSplineIndexSet<GV> indexSet_;
};

}   // namespace Functions

}   // namespace Dune

#endif   // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH