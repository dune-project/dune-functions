#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>
#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>

namespace Dune {
  namespace Functions {


    namespace {


      template<typename size_type>
      struct ClearSizeVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void pre(Node& node, TreePath treePath)
        {
          leaf(node,treePath);
          node.setSize(0);
        }

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          node.setOffset(offset_);
        }

        ClearSizeVisitor(size_type offset)
          : offset_(offset)
        {}

        const size_type offset_;

      };


      template<typename Entity, typename size_type>
      struct BindVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void pre(Node& node, TreePath treePath)
        {
          node.setOffset(offset_);
        }

        template<typename Node, typename TreePath>
        void post(Node& node, TreePath treePath)
        {
          node.setSize(offset_ - node.offset());
        }

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          node.setOffset(offset_);
          node.bind(entity_);
          offset_ += node.size();
        }

        BindVisitor(const Entity& entity, size_type offset = 0)
          : entity_(entity)
          , offset_(offset)
        {}

        const Entity& entity_;
        size_type offset_;

      };


      template<typename size_type>
      struct InitializeTreeVisitor :
        public TypeTree::TreeVisitor,
        public TypeTree::DynamicTraversal
      {
        template<typename Node, typename TreePath>
        void pre(Node& node, TreePath treePath)
        {
          node.setTreeIndex(treeIndex_);
          ++treeIndex_;
        }

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          node.setTreeIndex(treeIndex_);
          ++treeIndex_;
        }

        InitializeTreeVisitor(size_type treeIndexOffset = 0) :
          treeIndex_(treeIndexOffset)
        {}

        size_type treeIndex_;
      };

    }


    template<typename TP>
    class BasisNodeMixin
    {

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,typename>
      friend struct BindVisitor;

      template<typename>
      friend struct InitializeTreeVisitor;

    public:

      using TreePath = TP;
      using size_type = std::size_t;

      BasisNodeMixin(const TreePath& treePath) :
        offset_(0),
        size_(0),
        treePath_(treePath),
        treeIndex_(0)
      {}

      size_type localIndex(size_type i) const
      {
        return offset_ + i;
      }

      size_type size() const
      {
        return size_;
      }

      const TreePath& treePath() const
      {
        return treePath_;
      }

      const size_type treeIndex() const
      {
        return treeIndex_;
      }

      size_type offset() const
      {
        return offset_;
      }

    protected:

      void setOffset(const size_type offset)
      {
        offset_ = offset;
      }

      void setSize(const size_type size)
      {
        size_ = size;
      }

      void setTreeIndex(const size_type treeIndex)
      {
        treeIndex_ = treeIndex;
      }

    private:

      size_type offset_;
      size_type size_;
      const TreePath treePath_;
      size_type treeIndex_;

    };


    template<typename SIZE_T_DUMMY, typename TP>
    class LeafBasisNode :
        public BasisNodeMixin<TP>,
        public TypeTree::LeafNode
    {

      using Mixin = BasisNodeMixin<TP>;

    public:

      using TreePath = TP;
      using size_type = std::size_t;

      LeafBasisNode(TreePath treePath = TreePath()) :
        Mixin(treePath)
      {}

    };


    template<typename SIZE_T_DUMMY, typename TP, typename T, std::size_t n>
    class PowerBasisNode :
      public BasisNodeMixin<TP>,
      public TypeTree::PowerNode<T,n>
    {

      using Mixin = BasisNodeMixin<TP>;
      using Node = TypeTree::PowerNode<T,n>;

    public:

      PowerBasisNode(const TP& tp) :
        Mixin(tp)
      {}

      PowerBasisNode(const TP& tp, const typename Node::NodeStorage& children) :
        Mixin(tp),
        Node(children)
      {}

    };


    template<typename SIZE_T_DUMMY, typename TP, typename... T>
    class CompositeBasisNode :
      public BasisNodeMixin<TP>,
      public TypeTree::CompositeNode<T...>
    {

      using Mixin = BasisNodeMixin<TP>;
      using Node = TypeTree::CompositeNode<T...>;

    public:

      CompositeBasisNode(const TP& tp)
        : Mixin(tp)
      {}

      CompositeBasisNode(const TP& tp, const typename Node::NodeStorage& children) :
        Mixin(tp),
        Node(children)
      {}

      template<typename... Children>
      CompositeBasisNode(const shared_ptr<Children>&... children, const TP& tp)
        : Mixin(tp)
        , Node(children...)
      {}

    };


    template<typename Tree>
    void clearSize(Tree& tree, std::size_t offset)
    {
      TypeTree::applyToTree(tree,ClearSizeVisitor<std::size_t>(offset));
    }

    template<typename Tree, typename Entity>
    void bindTree(Tree& tree, const Entity& entity, std::size_t offset = 0)
    {
      BindVisitor<Entity,std::size_t> visitor(entity,offset);
      TypeTree::applyToTree(tree,visitor);
    }

    template<typename Tree>
    void initializeTree(Tree& tree, std::size_t treeIndexOffset = 0)
    {
      InitializeTreeVisitor<std::size_t> visitor(treeIndexOffset);
      TypeTree::applyToTree(tree,visitor);
    }


  } // namespace Functions

} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
