// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH

#include <cassert>
#include <memory>
#include <vector>
#include <array>
#include <optional>

#include <dune/common/indices.hh>
#include <dune/common/tuplevector.hh>
#include <dune/common/typelist.hh>

#include <dune/typetree/hybridmultiindex.hh>
#include <dune/typetree/nodetags.hh>
#include <dune/typetree/traversal.hh>

namespace Dune {
  namespace Functions {


    namespace Impl {

      // This class encapsulates the access to the setOffset()
      // and setTreeIndex() methods of a node. This way we
      // can hide the methods from the user but still provide
      // access where this is needed.
      struct BasisNodeSetupHelper
      {

        template<class Node, class size_type>
        static void setSize(Node& node, const size_type size)
        {
          node.setSize(size);
        }

        template<class Node, class size_type>
        static void setOffset(Node& node, const size_type offset)
        {
          node.setOffset(offset);
        }

        template<class Node, class size_type>
        static void setTreeIndex(Node& node, const size_type index)
        {
          node.setTreeIndex(index);
        }

      };



      // A mixin class for generalized child access from
      // multiple indices or a tree path. The derived class
      // only has to provide the child(i) method with
      // a single index for accessing direct children.
      template<class Impl>
      class ChildAccessMixIn
      {

        Impl& asImpl()
        {
          return static_cast<Impl&>(*this);
        }

        const Impl& asImpl() const
        {
          return static_cast<const Impl&>(*this);
        }

      public:

        /**
         * \brief Const access to descendent node by indices
         *
         * \param ii Indices of descendents
         */
        template<class... II>
        const auto& child(II... ii) const
        requires (sizeof...(II) != 1)
        {
          return Dune::TypeTree::child(asImpl(), ii...);
        }

        /**
         * \brief Mutable access to descendent node by indices
         *
         * \param ii Indices of descendents
         */

        template<class... II>
        auto& child(II... ii)
        requires (sizeof...(II) != 1)
        {
          return Dune::TypeTree::child(asImpl(), ii...);
        }

        /**
         * \brief Const access to descendent node by tree path
         *
         * \param treePath Tree path identifying the descendent
         */
        template<class... II>
        const auto& child(Dune::HybridMultiIndex<II...> treePath) const
        {
          return Dune::TypeTree::child(asImpl(), treePath);
        }

        /**
         * \brief Mutable access to descendent node by tree path
         *
         * \param treePath Tree path identifying the descendent
         */
        template<class... II>
        auto& child(Dune::HybridMultiIndex<II...> treePath)
        {
          return Dune::TypeTree::child(asImpl(), treePath);
        }

      };

    } // end namespace Impl


    class BasisNodeMixin
    {

      friend struct Impl::BasisNodeSetupHelper;

    public:

      using size_type = std::size_t;

      BasisNodeMixin() :
        offset_(0),
        size_(0),
        treeIndex_(0)
      {}

      size_type localIndex(size_type i) const
      {
        assert(i < size_);
        return offset_ + i;
      }

      /**
       * \brief Obtain the number of basis function in the local node.
       *
       * Notice that it is undefined behaviour to access the `element()`
       * and `finiteElement()` methods of the node if it is empty, i.e.,
       * if its size is zero.
       */
      size_type size() const
      {
        return size_;
      }

      /**
       * \brief Check if the node is empty
       *
       * This is equivalent to `size()==0`.
       * Notice that it is undefined behaviour to access the `element()`
       * and `finiteElement()` methods of the node if it is empty, i.e.,
       * if its size is zero.
       */
      bool empty() const
      {
        return (size_ == 0);
      }

      size_type treeIndex() const
      {
        return treeIndex_;
      }

    protected:

      size_type offset() const
      {
        return offset_;
      }

      void setOffset(const size_type offset)
      {
        offset_ = offset;
      }

      void setSize(const size_type size)
      {
        size_ = size;
      }

      void setTreeIndex(size_type treeIndex)
      {
        treeIndex_ = treeIndex;
      }

    private:

      size_type offset_;
      size_type size_;
      size_type treeIndex_;

    };


    class LeafBasisNode :
        public BasisNodeMixin
    {
    public:

      // Begin of node interface

      static constexpr auto degree()
      {
        return Dune::index_constant<0>{};
      }

      // Historic node interface

      static const bool isLeaf [[deprecated]] = true;
      static const bool isPower [[deprecated]] = false;
      static const bool isComposite [[deprecated]] = false;
      using NodeTag [[deprecated]] = Dune::TypeTree::LeafNodeTag;

      // End of node interface

    };



    template<typename Node, typename Element>
    class InnerBasisNodeMixin
      : public BasisNodeMixin
    {
    public:

      void bind(const Element& entity)
      {
        Node& self = *static_cast<Node*>(this);
        std::size_t offset = this->offset();
        Dune::Hybrid::forEach(Dune::range(self.degree()), [&](auto i) {
          bindTree(self.child(i), entity, offset);
          offset += self.child(i).size();
        });
        this->setSize(offset - this->offset());
      }

    };



    template<typename T, std::size_t n>
    class PowerBasisNode :
      public InnerBasisNodeMixin<PowerBasisNode<T, n>, typename T::Element>,
      public Impl::ChildAccessMixIn<PowerBasisNode<T, n>>
    {
    public:

      // Begin of node interface

      static constexpr auto degree()
      {
        return Dune::index_constant<n>{};
      }

      template<class Index>
      requires (std::is_convertible_v<Index, std::size_t>)
      const auto& child(Index i) const
      {
        return children_[i].value();
      }

      template<class Index>
      requires (std::is_convertible_v<Index, std::size_t>)
      auto& child(Index i)
      {
        return children_[i].value();
      }

      using Impl::ChildAccessMixIn<PowerBasisNode<T, n>>::child;

      // Historic node interface

      using ChildType = T;

      static const bool isLeaf [[deprecated]] = false;
      static const bool isPower [[deprecated]] = true;
      static const bool isComposite [[deprecated]] = false;
      using NodeTag [[deprecated]] = Dune::TypeTree::PowerNodeTag;

      // End of node interface

      using Element = typename T::Element;

      PowerBasisNode() = default;

      const Element& element() const
      {
        return child(Dune::Indices::_0).element();
      }

      template<class Index, class TT>
      void setChild(Index i, TT&& t)
      {
        children_[i].emplace(std::forward<TT>(t));
      }

    private:
      std::array<std::optional<T>, n> children_;
    };



    template<typename T>
    class DynamicPowerBasisNode :
      public InnerBasisNodeMixin<DynamicPowerBasisNode<T>, typename T::Element>,
      public Impl::ChildAccessMixIn<DynamicPowerBasisNode<T>>
    {
    public:

      // Begin of node interface

      std::size_t degree() const
      {
        return children_.size();
      }

      template<class Index>
      requires (std::is_convertible_v<Index, std::size_t>)
      const auto& child(Index i) const
      {
        return children_[i].value();
      }

      template<class Index>
      requires (std::is_convertible_v<Index, std::size_t>)
      auto& child(Index i)
      {
        return children_[i].value();
      }

      using Impl::ChildAccessMixIn<DynamicPowerBasisNode<T>>::child;

      // Historic node interface

      using ChildType = T;

      static const bool isLeaf [[deprecated]] = false;
      static const bool isPower [[deprecated]] = true;
      static const bool isComposite [[deprecated]] = false;
      using NodeTag [[deprecated]] = Dune::TypeTree::DynamicPowerNodeTag;

      // End of node interface

      using Element = typename T::Element;

      DynamicPowerBasisNode (std::size_t children)
        : children_(children)
      {}

      const Element& element() const
      {
        return child(Dune::Indices::_0).element();
      }

      template<class Index, class TT>
      void setChild(Index i, TT&& t)
      {
        children_[i].emplace(std::forward<TT>(t));
      }

    private:
      std::vector<std::optional<T>> children_;
    };


    template<typename... T>
    class CompositeBasisNode :
      public InnerBasisNodeMixin<CompositeBasisNode<T...>, typename TypeListEntry_t<0, TypeList<T...>>::Element>,
      public Impl::ChildAccessMixIn<CompositeBasisNode<T...>>
    {
    public:

      // Begin of node interface

      static constexpr auto degree()
      {
        return Dune::index_constant<sizeof...(T)>{};
      }

      template<std::size_t i>
      const auto& child(Dune::index_constant<i> ii) const
      {
        return children_[ii].value();
      }

      template<std::size_t i>
      auto& child(Dune::index_constant<i> ii)
      {
        return children_[ii].value();
      }

      using Impl::ChildAccessMixIn<CompositeBasisNode<T...>>::child;

      // Historic node interface

      using ChildTypes = std::tuple<T...>;

      static const bool isLeaf [[deprecated]] = false;
      static const bool isPower [[deprecated]] = false;
      static const bool isComposite [[deprecated]] = true;
      using NodeTag [[deprecated]] = Dune::TypeTree::CompositeNodeTag;

      template<std::size_t k>
      struct Child {
        static_assert((k < degree()), "child index out of range");

        //! The type of the child.
        using Type = typename std::tuple_element_t<k,ChildTypes>;

        using type = Type;
      };

      // End of node interface

      using Element = typename Child<0>::Type::Element;

      CompositeBasisNode() = default;

      explicit CompositeBasisNode(const T&... children) :
        children_(children...)
      {}

      const Element& element() const
      {
        return child(Dune::Indices::_0).element();
      }

      template<std::size_t i, class TT>
      void setChild (TT&& t, Dune::index_constant<i> ii = {})
      {
        children_[ii].emplace(std::forward<TT>(t));
      }

    private:
      Dune::TupleVector<std::optional<T>...> children_;
    };


    template<typename Tree, typename Entity>
    void bindTree(Tree& tree, const Entity& entity, std::size_t offset = 0)
    {
      Impl::BasisNodeSetupHelper::setOffset(tree, offset);
      tree.bind(entity);
    }

    template<typename Tree>
    void initializeTree(Tree& tree, std::size_t treeIndexOffset = 0)
    {
      Dune::TypeTree::forEachNode(tree, [&](auto& node, const auto& treePath) {
          Impl::BasisNodeSetupHelper::setTreeIndex(node, treeIndexOffset);
          ++treeIndexOffset;
        });
    }


  } // namespace Functions

} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
