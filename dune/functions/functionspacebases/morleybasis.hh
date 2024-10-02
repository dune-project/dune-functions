// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MORLEYBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MORLEYBASIS_HH

#include <algorithm>
#include <type_traits>
#include <vector>
#include <array>
#include <bitset>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/functions/analyticfunctions/monomialset.hh>
#include <dune/functions/functionspacebases/cubichermitebasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/functionaldescriptor.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/transformedfiniteelementmixin.hh>

/**
 * \file morleybasis.hh
 * \brief This file provides an implementation of the quadratic Morley finite element in 1 to 3 dimensions.
 *
 * The implementation is based on [Robert Kirby, A general approach to transforming
 * finite elements, 2018].
 * It contains in the following order:
 *     - A GlobalBasis typedef MorleyBasis
 *     - A template H2LocalBasisTraits, extending the dune-localfunctions
 *       LocalBasisTraits by an exported Hessiantype
 *     - A template MorleyLocalFiniteElement providing an implementation
 *       of the LocalFiniteElement interface, along with its subparts (Impl namespace)
 *     - A template MorleyNode
 *     - A template MorleyPreBasis
 *     - A factory morley() in the BasisFactory namespace
 */
namespace Dune::Functions
{

  template<class GV, class R>
  struct MorleyPreBasis;

  /** \brief Nodal basis of a scalar quadratic Morley finite element space
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \note This only works for simplex grids. The Morley basis is only implemented for degree 2.
   * \note The Morley Finite element has the following properties:
   *   - Its global space is not in \f$ H^1 \f$ but can be used for nonconforming method in \f$ H^2\f$.
   *   - Its interpolation evaluates derivatives, i.e. you cannot interpolate into a lambda function.
   *   - Strongly enforcing boundary conditions is not as simple as with langrange bases.
   *   - It global space is not nested, i.e. the space on a refined grid is not a subspace of the
   *     space on the coarser grid.
   * All arguments passed to the constructor will be forwarded to the constructor
   * of MorleyPreBasis.
   *
   * \tparam GV The GridView that the space is defined on
   * \tparam R The range type of the local basis
   */
  template<class GV, class R = double>
  using MorleyBasis = DefaultGlobalBasis<MorleyPreBasis<GV, R> >;


  namespace Impl
  {

    /**
     * \brief Implementation of morley Polynomials
     * \tparam D Type to represent the field in the domain
     * \tparam R Type to represent the field in the range
     * \tparam dim Dimension of the domain simplex
     */
    template<class D, class R>
    class MorleyReferenceLocalBasis
    {
    public:
      static constexpr int dim = 2;
      using Traits = H2LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>,
                                        FieldMatrix<R, 1, dim>, FieldMatrix<R, dim, dim>>;

    private:

      /**
       * \brief Get the Morley Coefficients Matrix
       * \return FieldMatrix<F, size, size>
       *  where size is the dimension of the cubic polynomial space
       *
       * This returns the basis transformation matrix from a monomial
       * basis to the Morley basis on the reference domain.
       * I.e. the i-th row of the returned matrix contains the
       * coefficients of the i-th Morley basis function with respect to a
       * monomial basis. The monomials are enumerated as in the
       * MonomialSet of the corresponding dimension and order.
       */
      static constexpr auto getMorleyCoefficients()
      {
        // Define std::sqrt(2.) manually in double precision,
        // because std::sqrt is not constexpr before C++26.
        D sqrt2 = 0.5 * 1.414213562373095;

        return Dune::FieldMatrix<D, 6,6>({{1, -1, -1, 0, 2, 0},
                                {0, 0.5, 0.5, 0.5, -1, -0.5},
                                {0, 0.5, 0.5, -0.5, -1, 0.5},
                                {0, 0, 1, 0, 0, -1},
                                {0, -1., 0, 1., 0, 0},
                                {0, sqrt2, sqrt2, -sqrt2, -2. * sqrt2, -sqrt2}});
      }

      static constexpr auto referenceBasisCoefficients = getMorleyCoefficients();
      MonomialSet<typename Traits::RangeFieldType, dim, 2> monomials;

    public:

      /** The number of basis functions in the basis
       */
      static constexpr unsigned int size()
      {
        return 6;
      }

      /** The polynomial order of the basis
       */
      static constexpr unsigned int order()
      {
        return 2;
      }

      /** \brief Evaluate function values of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Values of all shape functions at that point
       */
      void evaluateFunction(const typename Traits::DomainType &in,
                            std::vector<typename Traits::RangeType> &out) const
      {
        out.resize(size());
        auto monomialValues = monomials(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
      }

      /** \brief Evaluate Jacobians of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Jacobians of all shape functions at that point
       */
      void evaluateJacobian(const typename Traits::DomainType &in,
                            std::vector<typename Traits::JacobianType> &out) const
      {
        out.resize(size());
        auto monomialValues = derivative(monomials)(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
      }

      /** \brief Evaluate Hessians of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Hessians of all shape functions at that point
       */
      void evaluateHessian(const typename Traits::DomainType &in,
                            std::vector<typename Traits::HessianType> &out) const
      {
        out.resize(size());
        auto monomialValues = derivative(derivative(monomials))(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
      }

      /** \brief Evaluate partial derivatives of all shape functions at a given point
       *
       * \param[in] order The partial derivative to be computed, as a multi-index
       * \param[in] in  The evaluation point
       * \param[out] out Jacobians of all shape functions at that point
       */
      void partial(std::array<unsigned int, dim> order, const typename Traits::DomainType &in,
                  std::vector<typename Traits::RangeType> &out) const
      {
        out.resize(size());
        auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
        if (totalOrder == 0)
          evaluateFunction(in, out);
        else if (totalOrder == 1)
        {
          evaluateJacobian(in,jacobiansBuffer_);
          std::size_t which = std::max_element(order.begin(), order.end()) - order.begin();
          for (auto i : Dune::range(size()))
            out[i] = jacobiansBuffer_[i][0][which];
        }
        else if (totalOrder == 2)
        {
          evaluateHessian(in, hessianBuffer_);
          std::size_t first, second;
          first = std::max_element(order.begin(), order.end()) - order.begin();
          if (order[first] == 2)
            second = first;
          else
          {
            order[first] = 0;
            second = std::max_element(order.begin(), order.end()) - order.begin();
          }
          for (auto i : Dune::range(size()))
            out[i] = hessianBuffer_[i][first][second];
        }
        else
          DUNE_THROW(RangeError, "partial() not implemented for given order");
      }

    private:
      mutable std::vector<typename Traits::JacobianType> jacobiansBuffer_;
      mutable std::vector<typename Traits::HessianType> hessianBuffer_;
    };


    /** \brief Associations of the Morley degrees of freedom to subentities of the
     * reference simplex
     */
    class MorleyLocalCoefficients
    {
    public:
      using size_type = std::size_t;

      MorleyLocalCoefficients()
      {
        for (size_type i = 0; i < 3; ++i)
        {
          localKeys_[i] = LocalKey(i, 2, 0);     // vertex dofs
          localKeys_[3 + i] = LocalKey(i, 1, 0); // edge dofs
        }
      }

      /** number of coefficients
       */
      static constexpr size_type size()
      {
        return 6;
      }

      /** get i'th index
       */
      constexpr LocalKey const &localKey(size_type i) const
      {
        return localKeys_[i];
      }

    private:
      std::array<LocalKey, 6> localKeys_;
    };


    /**
     * \brief Class that evaluates the push forwards of the global nodes of a
     * LocalFunction. It stretches the LocalInterpolation interface, because we
     * evaluate the derivatives of f.
     */
    template<class D>
    class MorleyLocalInterpolation
    {
      using size_type = std::size_t;
      using FunctionalDescriptor = Dune::Functions::Impl::FunctionalDescriptor<2>;

    public:

      MorleyLocalInterpolation()
      {
        descriptors_[0] = FunctionalDescriptor();
        descriptors_[1] = FunctionalDescriptor();
        descriptors_[2] = FunctionalDescriptor();
        descriptors_[3] = FunctionalDescriptor(1);
        descriptors_[4] = FunctionalDescriptor(1);
        descriptors_[5] = FunctionalDescriptor(1);
      }

      /** \brief bind the Interpolation to an element and a localState.
       */
      template <class Element>
      void bind(const Element &element, const std::bitset<3> &edgeOrientation)
      {
        auto geometry = element.geometry();

        // get global Normals and midpoints
        auto refElement = Dune::referenceElement<double, 2>(geometry.type());
        for (std::size_t i = 0; i < 3; ++i)
        {
          localVertices_[i] = refElement.position(i, 2);

          localMidpoints_[i] = refElement.position(i, 1);
          std::size_t lower = (i == 2) ? 1 : 0;
          std::size_t upper = (i == 0) ? 1 : 2;

          auto edge = geometry.global(refElement.position(upper, 2))
                    - geometry.global(refElement.position(lower, 2));
          // normalize and orient
          edge /= edge.two_norm() * (edgeOrientation[i] ? -1. : 1.);
          // Rotation by pi/2. Note that Kirby rotates by 3*pi/2
          globalNormals_[i] = {-edge[1], edge[0]};
        }
      }

      /** \brief Evaluate a given function at the Lagrange nodes
       *
       * \tparam F Type of function to evaluate
       * \tparam C Type used for the values of the function
       * \param[in] f Function to evaluate
       * \param[out] out Array of function values
       */
      template <class F, class C>
      void interpolate(const F &f, std::vector<C> &out) const
      {
        out.resize(6);
        auto&& df = derivative(f);
        for (size_type i = 0; i < 3; ++i)
        {
          out[i] = f(localVertices_[i]);
          out[3 + i] = squeeze(df(localMidpoints_[i])).dot(globalNormals_[i]);
        }
      }

      const FunctionalDescriptor& functionalDescriptor(size_type i) const
      {
        return descriptors_[i];
      }

    protected:
      std::array<Dune::FieldVector<D, 2>, 3> globalNormals_;
      std::array<Dune::FieldVector<D, 2>, 3> localMidpoints_;
      std::array<Dune::FieldVector<D, 2>, 3> localVertices_;
      std::array<FunctionalDescriptor, 6> descriptors_;

    };

  } // namespace Impl

  template<class D, class R>
  using MorleyLocalBasisTraits = typename Impl::MorleyReferenceLocalBasis<D, R>::Traits;



  /** \brief Morley finite element for simplices, as defined on the reference Element.
   * Note, that this is a non affine-equivalent finite element, that requires an additional transformation to the relate reference basis with the pullbacks of global basis.
   * For more Details, see <dune/functions/functionspacebases/morleybasis.hh>.
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam dim dimension of the reference element
   */
  template<class D, class R>
  class MorleyLocalFiniteElement
    : public Impl::TransformedFiniteElementMixin<MorleyLocalFiniteElement<D,R>, MorleyLocalBasisTraits<D, R>>
  {
    using Base = Impl::TransformedFiniteElementMixin< MorleyLocalFiniteElement<D,R>, MorleyLocalBasisTraits<D, R>>;
    friend class Impl::TransformedLocalBasis<MorleyLocalFiniteElement<D,R>, MorleyLocalBasisTraits<D, R>>;
    static constexpr int dim = 2;
  public:

    /** \brief Export number types, dimensions, etc.
     */
    using size_type = std::size_t;
    using Traits = LocalFiniteElementTraits<
        Impl::TransformedLocalBasis<MorleyLocalFiniteElement<D,R>, MorleyLocalBasisTraits<D, R>>,
        Impl::MorleyLocalCoefficients,
        Impl::MorleyLocalInterpolation<D>>;

    /** \brief Returns the assignment of the degrees of freedom to the element
     * subentities
     */
    const typename Traits::LocalCoefficientsType &localCoefficients() const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType &localInterpolation() const
    {
      return interpolation_;
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type()
    {
      return GeometryTypes::simplex(dim);
    }

    /** The size of the transformed finite element.
     */
    static constexpr size_type size()
    {
      return 6;
    }

    /** Binds the Finite Element to an element.
      */
    template<class Element>
    void bind(std::bitset<3> const& data, Element const &e)
    {
      edgeOrientation_ = data;

      fillMatrix(e.geometry());
      interpolation_.bind(e, edgeOrientation_);
    }

  protected:

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    Impl::MorleyReferenceLocalBasis<D, R> const& referenceLocalBasis() const
    {
      return basis_;
    }

    /** Apply the transformation. Note that we do not distinguish for
     * Scalar/Vector/Matrix Type,
     * but only assume the Values to be Elements of a Vectorspace.
     * We assume random access containers.
     */
    template<class InputValues, class OutputValues>
    void transform(InputValues const& inValues, OutputValues& outValues) const
    {
      mat_.mv(inValues, outValues);
    }

  private:

    /**
     * \brief Fill the transformationmatrix \f$ M \f$
     *
     * \tparam Geometry  the Geometry class
     * \param geometry   the geometry of the element we are bound to
     */
    template<class Geometry>
    void fillMatrix(Geometry const &geometry)
    {
      std::array<R, 3> B_11;
      std::array<R, 3> B_12;
      std::array<R, 3> l_inv;

      std::array<Dune::FieldVector<R, 2>, 3> referenceTangents;
      std::array<Dune::FieldVector<R, 2>, 3> globalTangents;

      // By default, edges point from the vertex with the smaller index
      // to the vertex with the larger index.

      // get local and global Tangents
      auto refElement = Dune::referenceElement<double, 2>(geometry.type());
      auto x = refElement.position(0,0);
      for (std::size_t i = 0; i < 3; ++i)
      {
        std::size_t lower = (i == 2) ? 1 : 0;
        std::size_t upper = (i == 0) ? 1 : 2;
        auto edge = refElement.position(upper, 2) - refElement.position(lower, 2);

        referenceTangents[i] = edge / edge.two_norm();

        auto globalEdge = geometry.global(refElement.position(upper, 2))
                          - geometry.global(refElement.position(lower, 2));

        l_inv[i] = 1. / globalEdge.two_norm();
        globalTangents[i] = globalEdge * l_inv[i];
      }

      auto jacobianTransposed = geometry.jacobianTransposed(x);
      for (std::size_t i = 0; i < 3; ++i)
      {
        B_11[i] = -referenceTangents[i][1]
                    *(-globalTangents[i][1] * jacobianTransposed[0][0]
                      + globalTangents[i][0] * jacobianTransposed[0][1])
                  + referenceTangents[i][0]
                    *(-globalTangents[i][1] * jacobianTransposed[1][0]
                      + globalTangents[i][0] * jacobianTransposed[1][1]);
        B_12[i] = -referenceTangents[i][1]
                    *(globalTangents[i][0] * jacobianTransposed[0][0]
                      + globalTangents[i][1] * jacobianTransposed[0][1])
                  + referenceTangents[i][0]
                    *(globalTangents[i][0] * jacobianTransposed[1][0]
                      + globalTangents[i][1] * jacobianTransposed[1][1]);
      }

      // Actually setup matrix
      int sign = -1;
      for (std::size_t i = 0; i < 3; ++i)
      {
        mat_[i][i] = 1.;
        for (std::size_t j = 0; j < 3; ++j)
        {
          if (j != (2 - i)) // dune specific edge order
          {
            mat_[j][3 + i] = sign * B_12[i] * l_inv[i];
            sign *= -1;
          }
        }
        mat_[3 + i][3 + i] = (edgeOrientation_[i] ? -1. : 1.) * B_11[i];
      }
    }

    // a finite element consists of a basis, coeffiecents and an interpolation
    typename Impl::MorleyReferenceLocalBasis<D, R> basis_;
    typename Traits::LocalCoefficientsType coefficients_;
    typename Traits::LocalInterpolationType interpolation_;
    // This is the matrix M in Kirbys paper
    Dune::FieldMatrix<R, 6, 6>mat_;
    // the local state, i.e. a collection of global information restricted to this element
    std::bitset<3> edgeOrientation_;

  };

  // *****************************************************************************
  // This is the reusable part of the basis. It contains
  //
  //   MorleyPreBasis
  //   MorleyNode
  //
  // The pre-basis allows to create the others and is the owner of possible shared
  // state. These components do _not_ depend on the global basis and local view
  // and can be used without a global basis.
  // *****************************************************************************

  template<class GV, class R>
  class MorleyNode
    : public LeafBasisNode
  {
    using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;
    using FiniteElement = MorleyLocalFiniteElement<typename GV::ctype, R>;


    MorleyNode(Mapper const& m, std::vector<std::bitset<3>> const& data)
      : mapper_(&m), data_(&data)
    {
      this->setSize(finiteElement_.size());
    }

    //! Return current element, throw if unbound
    Element const &element() const
    {
      return *element_;
    }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the
     * dune-localfunctions module
     */
    FiniteElement const &finiteElement() const
    {
      return finiteElement_;
    }

    //! Bind to element.
    void bind(Element const &e)
    {
      element_ = &e;
      finiteElement_.bind((*data_)[mapper_->index(e)], *element_);
    }

    //! The order of the local basis.
    unsigned int order() const
    {
      return finiteElement_.localBasis().order();
    }

  protected:
    FiniteElement finiteElement_;
    Element const* element_;
    Mapper const* mapper_;
    std::vector<std::bitset<3>> const* data_;
  };


  /**
   * \brief A pre-basis for a Morleybasis
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV  The grid view that the FE basis is defined on
   * \tparam R   Range type used for shape function values
   * \note This only works for simplex grids
   */
  template<class GV, class R>
  class MorleyPreBasis
    : public LeafPreBasisMapperMixin<GV>
  {
    using Base = LeafPreBasisMapperMixin<GV>;
    using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
    using Element = typename GV::template Codim<0>::Entity;
    using D = typename GV::ctype;
    static const std::size_t dim = GV::dimension;

    // helper methods to assign each subentity the number of dofs. Used by the LeafPreBasisMapperMixin.
    static constexpr auto morleyMapperLayout(Dune::GeometryType type, int gridDim)
    {
      assert(gridDim == 2);
      if (type.isVertex())
        return 1; // one evaluation dof per vertex
      if (type.isLine())
        return 1;
      if ((type.isTriangle()) )
        return 0;
      else
        return 0;
    }

  public:
    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

    //! Template mapping root tree path to type of created tree node
    using Node = MorleyNode<GridView, R>;

    //! Constructor for a given grid view object
    MorleyPreBasis(const GV &gv)
      : Base(gv, morleyMapperLayout)
      , mapper_({gv, mcmgElementLayout()})
    {
      updateState(gv);
    }

    //! Update the stored grid view, to be called if the grid has changed
    void update(GridView const &gv)
    {
      Base::update(gv);
      updateState(gv);
    }

    /**
     * \brief Create tree node
     */
    Node makeNode() const
    {
      return Node{mapper_, data_};
    }

  protected:

    void updateState(GridView const &gridView)
    {
      data_.resize(gridView.size(0));
      // compute orientation for all elements
      unsigned short orientation = 0;
      auto const& idSet = gridView.grid().globalIdSet();

      for (const auto &element : elements(gridView))
      {
        const auto &refElement = referenceElement(element);
        auto elementIndex = mapper_.index(element);

        orientation = 0;

        for (std::size_t i = 0; i < element.subEntities(dim - 1); i++)
        {
          // Local vertex indices within the element
          auto localV0 = refElement.subEntity(i, dim - 1, 0, dim);
          auto localV1 = refElement.subEntity(i, dim - 1, 1, dim);

          // Global vertex indices within the grid
          auto globalV0 = idSet.subId(element, localV0, dim);
          auto globalV1 = idSet.subId(element, localV1, dim);

          // The edge is flipped if the local ordering disagrees with global ordering
          if ((localV0 < localV1 && globalV0 > globalV1)
              || (localV0 > localV1 && globalV0 < globalV1))
          {
            orientation |= (1 << i);
          }
        }
        data_[elementIndex] = orientation;
      }
    }

    SubEntityMapper mapper_;
    std::vector<std::bitset<3>> data_;

  }; // class MorleyPreBasis

  namespace BasisFactory
  {

    /**
    * \brief construct a PreBasisFactory for the Morley Finite Element
    *
    * \tparam R RangeFieldType
    * \return the PreBasisFactory
    */
    template<class R = double>
    auto morley()
    {
      return [=](auto const &gridView) {
        return MorleyPreBasis<std::decay_t<decltype(gridView)>, R>(gridView);
      };
    }

  } // end namespace BasisFactory

} // end namespace Dune::Functions

#endif
