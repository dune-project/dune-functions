// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH

#include <cassert>
#include <memory>

#include <dune/common/indices.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/compositenode.hh>
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
        public BasisNodeMixin,
        public TypeTree::LeafNode
    {};



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
      public TypeTree::PowerNode<T,n>
    {

      using Node = TypeTree::PowerNode<T,n>;

    public:

      using Element = typename T::Element;

      PowerBasisNode() = default;

      PowerBasisNode(const typename Node::NodeStorage& children) :
        Node(children)
      {}

      const Element& element() const
      {
        return this->child(Dune::Indices::_0).element();
      }

    };



    template<typename T>
    class DynamicPowerBasisNode :
      public InnerBasisNodeMixin<DynamicPowerBasisNode<T>, typename T::Element>,
      public TypeTree::DynamicPowerNode<T>
    {

      using Node = TypeTree::DynamicPowerNode<T>;

    public:

      using Element = typename T::Element;

      DynamicPowerBasisNode (std::size_t children)
        : Node(children)
      {}

      DynamicPowerBasisNode (typename Node::NodeStorage children)
        : Node(std::move(children))
      {}

      const Element& element() const
      {
        return this->child(0).element();
      }

    };


    template<typename... T>
    class CompositeBasisNode :
      public InnerBasisNodeMixin<CompositeBasisNode<T...>, typename TypeTree::CompositeNode<T...>::template Child<0>::Type::Element>,
      public TypeTree::CompositeNode<T...>
    {

      using Node = TypeTree::CompositeNode<T...>;

    public:

      using Element = typename Node::template Child<0>::Type::Element;

      CompositeBasisNode() = default;

      CompositeBasisNode(const typename Node::NodeStorage& children) :
        Node(children)
      {}

      explicit CompositeBasisNode(const T&... children) :
        Node(children...)
      {}

      template<typename... Children>
      explicit CompositeBasisNode(const std::shared_ptr<Children>&... children) :
        Node(children...)
      {}

      const Element& element() const
      {
        return this->child(Dune::Indices::_0).element();
      }

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
