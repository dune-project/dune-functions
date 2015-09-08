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

        ComputeSizeVisitor(const Entity& entity, size_type offset = 0)
          : entity_(entity)
          , offset_(offset)
        {}

        const Entity& entity_;
        size_type offset_;

      };

    }


    template<typename size_t, typename TP>
    class LeafBasisNode
      : public TypeTree::LeafNode
    {

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,typename>
      friend struct BindVisitor;

    public:

      using TreePath = TP;
      using size_type = size_t;

      LeafBasisNode(TreePath treePath = TreePath()) :
        offset_(0),
        treePath_(treePath)
      {}

      size_type localIndex(size_type i) const
      {
        return offset_ + i;
      }

      const TreePath& treePath() const
      {
        return treePath_;
      }

    protected:

      void setOffset(const size_type offset)
      {
        offset_ = offset;
      }

    private:

      size_type offset_;
      const TreePath treePath_;

    };


    template<typename size_t, typename TP>
    class InternalBasisNodeMixin
    {

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,typename>
      friend struct BindVisitor;

    public:

      using TreePath = TP;
      using size_type = size_t;

      InternalBasisNodeMixin() :
        offset_(0),
        size_(0)
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

    protected:

      void setOffset(const size_type offset)
      {
        offset_ = offset;
      }

      void setSize(const size_type size)
      {
        size_ = size;
      }

    private:

      size_type offset_;
      size_type size_;

    };


    template<typename size_t, typename TP, typename T, std::size_t n>
    class PowerBasisNode :
      public InternalBasisNodeMixin<size_t,TP>,
      public TypeTree::PowerNode<T,n>
    {

      using Mixin = InternalBasisNodeMixin<size_t,TP>;
      using Node = TypeTree::PowerNode<T,n>;

    protected:

      PowerBasisNode(const TP& tp) :
        Mixin(tp)
      {}

      PowerBasisNode(const TP& tp, const typename Node::NodeStorage& children) :
        Mixin(tp),
        Node(children)
      {}

    };


    template<typename size_t, typename TP, typename... T>
    class CompositeBasisNode :
      public InternalBasisNodeMixin<size_t,TP>,
      public TypeTree::CompositeNode<T...>
    {

      using Mixin = InternalBasisNodeMixin<size_t,TP>;
      using Node = TypeTree::CompositeNode<T...>;

    protected:

      CompositeBasisNode(const TP& tp)
        : Mixin(tp)
      {}

      CompositeBasisNode(const TP& tp, const typename Node::NodeStorage& children) :
        Mixin(tp),
        Node(children)
      {}

      template<typename... Children>
      CompositeBasisNode(const shared_ptr<Children>&... children)
        : Mixin(tp)
        , Node(children...)
      {}

    };


    template<typename Tree, typename size_type>
    void clearSize(Tree& tree, size_type offset)
    {
      TypeTree::applyToTree(tree,ClearSizeVisitor<size_type>(offset));
    }

    template<typename Tree, typename Entity, typename size_type = std::size_t>
    void bindTree(Tree& tree, const Entity& entity, size_type offset = 0)
    {
      BindVisitor<Entity,size_type> visitor(entity,offset);
      TypeTree::applyToTree(tree,visitor);
    }


  } // namespace Functions

} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
