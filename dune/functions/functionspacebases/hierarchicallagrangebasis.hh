// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICALLAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICALLAGRANGEBASIS_HH

#include <dune/common/exceptions.hh>
#include <dune/localfunctions/hierarchical/hierarchicalp2.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {
  namespace Functions {

    // *****************************************************************************
    // Implementation for Hierarchical Lagrange Basis
    //
    // -- only order k=2 is implemented up to now --
    // -- currently only supports simplex grids --
    //
    // This is the reusable part of the HierarchicalLagrangeBasis. It contains
    //
    //  HierarchicalLagrangePreBasis
    //  HierarchicalLagrangeNode
    //
    // The pre-basis allows to create the others and is the owner of possible shared
    // state. These components do _not_ depend on the global basis and can be
    // used without a global basis.
    // *****************************************************************************

    template<typename GV, int k, typename R=double>
    class HierarchicalLagrangeNode;

    template<typename GV, int k, class MI, typename R=double>
    class HierarchicalLagrangePreBasis;

    /**
     * \brief A pre-basis for a hierarchical basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV  The grid view that the FE basis is defined on
     * \tparam k   The polynomial order of ansatz functions (limited to second order till now)
     * \tparam MI  Type to be used for multi-indices
     * \tparam R   Range type used for shape function values
     */
    template<typename GV, int k, class MI, typename R>
    class HierarchicalLagrangePreBasis
    {
      static const int dim = GV::dimension;

    public:

      //! The grid view that the FE basis is defined on
      using GridView = GV;

      //! Type used for indices and size information
      using size_type = std::size_t;

      //! Template mapping root tree path to type of created tree node
      using Node = HierarchicalLagrangeNode<GV, k, R>;

      //! Type of created tree node index set
      using IndexSet = Impl::DefaultNodeIndexSet<HierarchicalLagrangePreBasis>;

      //! Type used for global numbering of the basis vectors
      using MultiIndex = MI;

      //! Type used for prefixes handed to the size() method
      using SizePrefix = Dune::ReservedVector<size_type, 1>;

      /** \brief Constructor for a given grid view object with layout for second order
       *
       *  (adjust for higher-orders)
       */
      HierarchicalLagrangePreBasis(const GridView& gv) : gridView_(gv) , mcmgMapper_(gv,p2Layout())
      {}

      //! Initialize the global indices
      void initializeIndices()
      {}

      //! Obtain the grid view that the basis is defined on
      const GridView& gridView() const
      {
        return gridView_;
      }

      //! Update the stored grid view & `MultipleCodimMultipleGeomTypeMapper`, to be called if the grid has changed
      void update (const GridView& gv)
      {
        gridView_ = gv;
        mcmgMapper_.update(gv);
      }

      /**
       * \brief Create tree node
       */
      Node makeNode() const
      {
        return Node{};
      }

      /**
       * \brief Create tree node index set
       *
       * Create an index set suitable for the tree node obtained
       * by makeNode().
       */
      IndexSet makeIndexSet() const
      {
        return IndexSet{*this};
      }

      //! Same as size(prefix) with empty prefix
      size_type size() const
      {
        return mcmgMapper_.size();
      }

      //! Return number of possible values for next position in multi index
      size_type size(const SizePrefix prefix) const
      {
        assert(prefix.size() == 0 || prefix.size() == 1);
        return (prefix.size() == 0) ? size() : 0;
      }

      //! Get the total dimension of the space spanned by this basis
      size_type dimension() const
      {
        return size();
      }

      /** \brief Get the maximal number of DOFs associated to node for any element
       *
       * See https://en.wikipedia.org/wiki/Figurate_number for an explanation of the formula
       */
      size_type maxNodeSize() const
      {
        // That cast to unsigned int is necessary because GV::dimension is an enum
        return Dune::binomial(std::size_t(order() + (unsigned int)GV::dimension),std::size_t(order()));
      }

      template<typename It>
      It indices(const Node& node, It it) const
      {
        for (size_type i = 0, end = node.finiteElement().size() ; i < end ; ++it, ++i)
        {
          Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
          const auto& element = node.element();

          *it = {{ (size_type)(mcmgMapper_.subIndex(element,localKey.subEntity(),localKey.codim())) }};
        }
        return it;
      }

    protected:
      GridView gridView_;

      unsigned int order() const
      {
        return 2;
      }

      MultipleCodimMultipleGeomTypeMapper<GridView> mcmgMapper_;

    private:
      /** \brief layout function for `MultipleCodimMultipleGeomTypeMapper`
       *
       *  (adjust for higher-orders)
       */
      static auto p2Layout()
      {
        return [](Dune::GeometryType type, int gridDim)
               {
                 if (type.isVertex())
                   return 1;
                 if (type.isLine())
                   return 1;
                 if (type.isTriangle())
                   return 0;
                 assert(type.isTetrahedron());
                 return 0;
               };
      }
    };



    template<typename GV, int k, typename R>
    class HierarchicalLagrangeNode :
      public LeafBasisNode
    {
      static const int dim = GV::dimension;

    public:

      using size_type = std::size_t;
      using Element = typename GV::template Codim<0>::Entity;
      using FiniteElement = HierarchicalP2LocalFiniteElement<typename GV::ctype,R,dim>;

      HierarchicalLagrangeNode() :
        finiteElement_(),
        element_(nullptr)
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

        if (e.type() != finiteElement_.type())
          DUNE_THROW(Dune::Exception,
                     "HierarchicalLagrange-elements do not exist for elements of type " << e.type());

        this->setSize(finiteElement_.size());
      }

    protected:

      unsigned int order() const
      {
        return 2;
      }

      const FiniteElement finiteElement_;
      const Element* element_;
    };



    namespace BasisFactory {

      namespace Impl {

        template<int k, typename R=double>
        class HierarchicalLagrangePreBasisFactory
        {

        public:
          static const std::size_t requiredMultiIndexSize = 1;

          template<class MultiIndex, class GridView>
          auto makePreBasis(const GridView& gridView) const
          {
            return HierarchicalLagrangePreBasis<GridView, k, MultiIndex, R>(gridView);
          }

        };
      } // end namespace BasisFactory::Impl


      /**
       * \brief Create a pre-basis factory that can create a HierarchicalLagrange pre-basis
       *
       * \ingroup FunctionSpaceBasesImplementations
       *
       * \tparam k   The polynomial order of the ansatz functions (limited to second order till now)
       * \tparam R   The range type of the local basis
       */
      template<std::size_t k, typename R=double>
      auto hierarchicalLagrange()
      {
        return Impl::HierarchicalLagrangePreBasisFactory<k,R>();
      }

    } // end namespace BasisFactory

    /** \brief Basis of a scalar Hierarchical Lagrange finite element space
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV The GridView that the space is defined on
     * \tparam k The order of the basis (limited to second order till now)
     * \tparam R The range type of the local basis
     *
     *  -- currently only supports simplex grids --
     */
    template<typename GV, int k, typename R=double>
    using HierarchicalLagrangeBasis = DefaultGlobalBasis<HierarchicalLagrangePreBasis<GV, k, FlatMultiIndex<std::size_t>, R> >;

  } // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICALLAGRANGEBASIS_HH
