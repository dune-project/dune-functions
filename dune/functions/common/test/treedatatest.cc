// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <memory>
#include <functional>

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/common/treedata.hh>

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

using namespace Dune;
using namespace Dune::Functions;


template <class Node>
using NodeData = double;

template <class Tree>
bool operator==(TreeData<Tree,NodeData,false> const& t1, TreeData<Tree,NodeData,false> const& t2)
{
  assert(t1.tree() == t2.tree() && t1.tree() != nullptr);

  bool same = true;
  forEachNode(*t1.tree(), [&](auto const& node, auto) {
    same = same && (t1[node] == t2[node]);
  });

  return same;
}

template <class Tree>
bool operator==(TreeData<Tree,NodeData,true> const& t1, TreeData<Tree,NodeData,true> const& t2)
{
  assert(t1.tree() == t2.tree() && t1.tree() != nullptr);

  bool same = true;
  forEachLeafNode(*t1.tree(), [&](auto const& node, auto) {
    same = same && (t1[node] == t2[node]);
  });

  return same;
}

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance(argc, argv);
  using namespace Functions::BasisBuilder;

  bool passed = true;

  const int dim = 2;
  typedef YaspGrid<2> GridType;
  FieldVector<double,2> bbox = {1, 1};
  std::array<int,2> num = {1, 1};
  GridType grid(bbox,num);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();


  auto taylorHoodBasis = makeBasis(
    gridView,
    composite(
      power<2>(lagrange<2>()),
      lagrange<1>()
    ));

  auto localView = taylorHoodBasis.localView();
  auto const& tree = localView.tree();
  using Tree = std::remove_reference_t<decltype(tree)>;

  // test treeData for all nodes
  {
    using TypeTree::forEachNode;

    // call default-constructor
    TreeData<Tree, NodeData, false> treeData;
    treeData.init(tree);

    forEachNode(tree, [&](auto const& node, auto) {
      treeData[node] = double(node.treeIndex());
    });

    // call constructor with tree
    TreeData<Tree, NodeData, false> treeData1(tree);

    // call init on non-empty treeData
    treeData1.init(tree);

    // call copy-constructor
    TreeData<Tree, NodeData, false> treeData2(treeData);
    passed = passed && (treeData == treeData2);

    // call copy-assignment operator on empty treeData
    TreeData<Tree, NodeData, false> treeData3;
    treeData3 = treeData;
    passed = passed && (treeData == treeData3);

    // call copy-assignment operator on non-empty treeData
    treeData2 = treeData3;
    passed = passed && (treeData3 == treeData2);

    // call move-assignment operator on non-empty treeData
    treeData = std::move(treeData2);

    // call move-constructor
    TreeData<Tree, NodeData, false> treeData4(std::move(treeData3));
  }

  // test treeData for leaf only
  {
    using TypeTree::forEachLeafNode;

    // call default-constructor
    TreeData<Tree, NodeData, true> treeData;
    treeData.init(tree);

    forEachLeafNode(tree, [&](auto const& node, auto) {
      treeData[node] = double(node.treeIndex());
    });

    // call constructor with tree
    TreeData<Tree, NodeData, true> treeData1(tree);

    // call init on non-empty treeData
    treeData1.init(tree);

    // call copy-constructor
    TreeData<Tree, NodeData, true> treeData2(treeData);
    passed = passed && (treeData == treeData2);

    // call copy-assignment operator on empty treeData
    TreeData<Tree, NodeData, true> treeData3;
    treeData3 = treeData;
    passed = passed && (treeData == treeData3);

    // call copy-assignment operator on non-empty treeData
    treeData2 = treeData3;
    passed = passed && (treeData3 == treeData2);

    // call move-assignment operator on non-empty treeData
    treeData = std::move(treeData2);

    // call move-constructor
    TreeData<Tree, NodeData, true> treeData4(std::move(treeData3));
  }

  // test for operations with uninitialized tree
  {
    // call default-constructor without initialization
    TreeData<Tree, NodeData, true> treeData;

    // call copy-constructor
    TreeData<Tree, NodeData, true> treeData2(treeData);

    // call move-constructor
    TreeData<Tree, NodeData, true> treeData3(std::move(treeData));
  }

  return passed ? 0: 1;
}
catch( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
