// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/grid/yaspgrid.hh>


namespace CD = Dune::Functions::ContainerDescriptors;

namespace TestImpl {

// check at run time whether index is a valid child index
template <class Node, class Index>
std::true_type checkAccessIndex (Node const& node, Index i)
{
  assert(std::size_t(i) < node.size() && "Child index out of range");
  return {};
}

// check at compile time whether index is a valid index
template <class Node, std::size_t i>
std::bool_constant<(i < Node::size())> checkAccessIndex (Node const& node, Dune::index_constant<i>)
{
  static_assert(i < Node::size(), "Child index out of range");
  return {};
}

// finally return the node itself if no further indices are provided. Break condition
// for the recursion over the node children.
template<class Node>
decltype(auto) accessImpl (Node&& node)
{
  return std::forward<Node>(node);
}

// recursively call `node[i]` with the given indices
template<class Node, class I0, class... I>
decltype(auto) accessImpl (Node&& node, I0 i0, [[maybe_unused]] I... i)
{
  auto valid = checkAccessIndex(node,i0);
  if constexpr (valid)
    return accessImpl(node[i0],i...);
  else
    return;
}

// forward to the impl methods by extracting the indices from the tree path
template<class Tree, class... Indices, std::size_t... i>
decltype(auto) access (Tree&& tree, [[maybe_unused]] Dune::TypeTree::HybridTreePath<Indices...> tp, std::index_sequence<i...>)
{
  return accessImpl(std::forward<Tree>(tree),Dune::TypeTree::treePathEntry<i>(tp)...);
}

// access a tree using a hybridTreePath
template<class Tree, class... Indices>
decltype(auto) access (Tree&& tree, Dune::TypeTree::HybridTreePath<Indices...> tp)
{
  return access(std::forward<Tree>(tree),tp,std::index_sequence_for<Indices...>{});
}

// convert a HybridTreePath into a ReservedVector
template <class SizePrefix, class... Indices>
auto sizePrefix (Dune::TypeTree::HybridTreePath<Indices...> tp)
{
  using value_type = typename SizePrefix::value_type;
  return Dune::unpackIntegerSequence([&](auto... i) {
    return SizePrefix{value_type(tp[i])...};
  }, std::index_sequence_for<Indices...>{});
}

} // end namespace TestImpl


// check that the sizes of an container descriptor correspond to the sizes provided by the
// basis (size-provider) directly.
template<class ContainerDescriptor, class SizeProvider,
         class PrefixPath = Dune::TypeTree::HybridTreePath<>>
void checkSize (Dune::TestSuite& test, const ContainerDescriptor& cd,
                const SizeProvider& sizeProvider, PrefixPath prefix = {})
{
  using SizePrefix = typename SizeProvider::SizePrefix;
  auto size1 = sizeProvider.size(TestImpl::sizePrefix<SizePrefix>(prefix));
  auto size2 = Dune::Hybrid::size(TestImpl::access(cd, prefix));
  test.require(std::size_t(size1) == std::size_t(size2), "size1 == size2");

  if constexpr(PrefixPath::size() < SizePrefix::max_size()) {
    Dune::Hybrid::forEach(Dune::range(size2), [&](auto i) {
      checkSize(test, cd, sizeProvider, push_back(prefix,i));
    });
  }
}

// check a specific multi-index by traversing all its components and the container descriptor
//  simultaneously
template <class ContainerDescriptor, class MultiIndex>
void checkMultiIndex (Dune::TestSuite& test, const ContainerDescriptor& cd,
                      const MultiIndex& mi, std::size_t j = 0)
{
  if (j < mi.size()) {
    test.check(mi[j] < cd.size(), "mi[j] < cd.size");
    auto size = Dune::Hybrid::size(cd);
    Dune::Hybrid::switchCases(Dune::range(size), mi[j],
      [&](auto jj) { checkMultiIndex(test,cd[jj],mi,j+1); });
  }
}

// check that all multi-indices of a global basis are within the range of the container descriptor
template <class ContainerDescriptor, class Basis>
void checkMultiIndices (Dune::TestSuite& test, const ContainerDescriptor& cd, const Basis& basis)
{
  auto localView = basis.localView();
  for (auto const& e : elements(basis.gridView()))
  {
    localView.bind(e);
    for (std::size_t i = 0; i < localView.size(); ++i) {
      auto mi = localView.index(i);
      checkMultiIndex(test,cd,mi);
    }
  }
}

void checkHierarchic (Dune::TestSuite& test)
{
  using namespace Dune::Indices;
  CD::Tuple<CD::Array<CD::FlatVector,3>,CD::FlatVector> stokes{
    CD::Array<CD::FlatVector,3>{CD::FlatVector{10},CD::FlatVector{10},CD::FlatVector{10}}, CD::FlatVector{5} };

  { // check the make-functions
    auto stokes2 = CD::makeDescriptor(
      CD::makeDescriptor(
        CD::makeUniformDescriptor(10, CD::Value{}),
        CD::makeUniformDescriptor(10, CD::Value{}),
        CD::makeUniformDescriptor(10, CD::Value{})
      ),
      CD::makeUniformDescriptor(5, CD::Value{})
    );

    static_assert(std::is_same_v<decltype(stokes),decltype(stokes2)>);
  }

  CD::Array<CD::FlatVector,3> const& velocity = stokes[_0];
  test.check(velocity[0].size() == 10, "v[0].size == 10");
  test.check(velocity[1].size() == 10, "v[1].size == 10");
  test.check(velocity[2].size() == 10, "v[2].size == 10");

  CD::FlatVector const& pressure = stokes[_1];
  test.check(pressure.size() == 5, "p.size == 5");
}

void checkConstructors (Dune::TestSuite& test)
{
  CD::Unknown cd0{};
  CD::Value cd1{};

  // the children might have different types
  CD::Tuple<CD::Unknown,CD::Value> cd2a;
  CD::Tuple<CD::Unknown,CD::Value> cd2b{cd0,cd1};
  auto cd2c = CD::makeDescriptor(cd0,cd1);
  static_assert(cd2a.size() == 2, "cd2a.size() == 2");
  static_assert(cd2b.size() == 2, "cd2b.size() == 2");
  static_assert(cd2c.size() == 2, "cd2c.size() == 2");
  static_assert(std::is_same_v<decltype(cd2c), decltype(cd2a)>);


  // all childresn have the same type, the number of children is static
  CD::Array<CD::Value, 3> cd3a;
  CD::Array<CD::Value, 3> cd3b = Dune::filledArray<3>(cd1);
  CD::Array<CD::Value, 3> cd3c{cd1,cd1,cd1};
  auto cd3d = CD::makeDescriptor(cd1,cd1,cd1);
  static_assert(cd3a.size() == 3, "cd3a.size() == 3");
  static_assert(cd3b.size() == 3, "cd3b.size() == 3");
  static_assert(cd3c.size() == 3, "cd3c.size() == 3");
  static_assert(cd3d.size() == 3, "cd3d.size() == 3");
  static_assert(std::is_same_v<decltype(cd3c), decltype(cd3a)>);
  static_assert(std::is_same_v<decltype(cd3d), decltype(cd3a)>);

  // all childresn have the same type, the number of children is a runtime value
  CD::Vector<CD::Value> cd4a(4,cd1);
  CD::Vector<CD::Value> cd4b{cd1,cd1,cd1,cd1};
  CD::Vector<CD::Value> cd4c(4,cd1);
  CD::Vector<decltype(cd1)> cd4d{cd1,cd1,cd1,cd1};
  test.check(cd4a.size() == 4, "cd4a.size() == 4");
  test.check(cd4b.size() == 4, "cd4b.size() == 4");
  test.check(cd4c.size() == 4, "cd4c.size() == 4");
  test.check(cd4d.size() == 4, "cd4d.size() == 4");
  static_assert(std::is_same_v<decltype(cd4c), decltype(cd4a)>);
  static_assert(std::is_same_v<decltype(cd4d), decltype(cd4a)>);

  // all children are identical, the number of children is static
  // only a single child is stored
  CD::UniformArray<CD::Value, 5> cd5a;
  CD::UniformArray<CD::Value, 5> cd5b(cd1);
  auto cd5c = makeUniformDescriptor(std::integral_constant<std::size_t,5>{},cd1);
  static_assert(cd5a.size() == 5, "cd5a.size() == 5");
  static_assert(cd5b.size() == 5, "cd5b.size() == 5");
  static_assert(cd5c.size() == 5, "cd5c.size() == 5");
  static_assert(std::is_same_v<decltype(cd5c), decltype(cd5a)>);

  // shortcut for uniform arrays storing `Value`
  CD::FlatArray<5> cd5e;
  static_assert(cd5e.size() == 5, "cd5e.size() == 5");
  static_assert(std::is_same_v<decltype(cd5e), decltype(cd5a)>);

  // all children are identical, the number of children is a runtime value
  // only a single child is stored
  CD::UniformVector<CD::Value> cd6a(6);
  CD::UniformVector cd6b(6,cd1);
  auto cd6c = makeUniformDescriptor(6,cd1);
  test.check(cd6a.size() == 6, "cd6a.size() == 6");
  test.check(cd6b.size() == 6, "cd6b.size() == 6");
  test.check(cd6c.size() == 6, "cd6c.size() == 6");
  static_assert(std::is_same_v<decltype(cd6b), decltype(cd6a)>);
  static_assert(std::is_same_v<decltype(cd6c), decltype(cd6a)>);

  // shortcut for uniform vectors storing `Value`
  CD::FlatVector cd6d(6);
  test.check(cd6d.size() == 6, "cd6d.size() == 6");
  static_assert(std::is_same_v<decltype(cd6d), decltype(cd6a)>);
}

template <class BasisFactory>
void checkBasis(Dune::TestSuite& test, const BasisFactory& bf)
{
  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0, 1.0}, {2, 2});
  auto basis = makeBasis(grid.leafGridView(), bf);
  checkSize(test, containerDescriptor(basis.preBasis()), basis);
  checkMultiIndices(test, containerDescriptor(basis.preBasis()), basis);
}

int main (int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;
  checkConstructors(test);
  checkHierarchic(test);


  using namespace Dune::Functions::BasisFactory;

  // test default index-merging strategies
  checkBasis(test, lagrange<1>() );
  checkBasis(test, power<2>(lagrange<1>()) );
  checkBasis(test, composite(lagrange<1>(),lagrange<2>()) );

  // test all combinations of two nested power bases
  Dune::Hybrid::forEach(std::tuple{flatLexicographic(), flatInterleaved(), blockedLexicographic(), blockedInterleaved()}, [&](auto outerIMS) {
    Dune::Hybrid::forEach(std::tuple{flatLexicographic(), flatInterleaved(), blockedLexicographic(), blockedInterleaved()}, [&](auto innerIMS) {
      checkBasis(test, power<2>(power<3>(lagrange<2>(), innerIMS), outerIMS) );
    });
  });

  // test more complicated bases
  checkBasis(test, power<2>(
        composite(
          power<1>(power<1>(lagrange<1>(), blockedInterleaved()), blockedLexicographic()),
          power<2>(lagrange<1>(), blockedInterleaved()),
          power<3>(lagrange<1>(), blockedLexicographic())
        ), blockedLexicographic()
      )
  );

  checkBasis(test, composite(
      composite(
        power<2>(lagrange<1>(), blockedInterleaved()),
        power<1>(power<2>(lagrange<1>(), blockedLexicographic()), blockedLexicographic())
      ),
      composite(
        power<1>(power<1>(lagrange<1>(), blockedLexicographic()), blockedInterleaved()),
        power<2>(lagrange<1>(), blockedInterleaved()),
        power<3>(lagrange<1>(), blockedLexicographic())
      )
    )
  );

  checkBasis(test, composite(
      composite(
        power<2>(lagrange<1>(), flatInterleaved()),
        power<1>(power<2>(lagrange<1>(), flatLexicographic()),flatLexicographic()),
        flatLexicographic()
      ),
      composite(
        power<1>(power<1>(lagrange<1>(), flatLexicographic()), flatInterleaved()),
        power<2>(lagrange<1>(), flatInterleaved()),
        power<3>(lagrange<1>(), flatLexicographic()),
        flatLexicographic()
      ),
      flatLexicographic()
    )
  );

  return test.exit();
}
