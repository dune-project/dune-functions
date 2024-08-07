// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_PYTHON_FUNCTIONS_TREE_HH
#define DUNE_PYTHON_FUNCTIONS_TREE_HH

#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>
#include <memory>

// for void_t, in <type_traits> for C++17
#include <dune/common/typetraits.hh>
#include <dune/common/visibility.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/python/localfunctions/localfiniteelement.hh>

namespace Dune {
namespace Python {
namespace Functions {

namespace detail {

template<typename Tree, std::enable_if_t< Tree::isComposite, int > = 0>
void registerTree_(pybind11::handle scope, const char* name = "Tree");

template<typename Tree, std::enable_if_t< Tree::isLeaf, int > = 0>
void registerTree_(pybind11::handle scope, const char* name = "Tree");

template<typename Tree, std::enable_if_t< Tree::isPower, int > = 0>
void registerTree_(pybind11::handle scope, const char* name = "Tree");

template<typename Tree>
void registerTreeCommon(pybind11::class_<Tree, std::shared_ptr<Tree>>& cls)
{
  /* dune-typetree properties */
  cls.def_property_readonly_static("isComposite", [](pybind11::object) { return Tree::isComposite; });
  cls.def_property_readonly_static("isLeaf", [](pybind11::object) { return Tree::isLeaf; });
  cls.def_property_readonly_static("isPower", [](pybind11::object) { return Tree::isPower; });
  cls.def("degree", [](pybind11::object) -> std::size_t { return Tree::degree(); });
  cls.def( "__len__", [] (const Tree& self) -> int { return self.size(); } );
  cls.def( "size", [] (const Tree& self) -> int { return self.size(); } );
  cls.def("localIndex", [](const Tree& self, unsigned int index) { return self.localIndex(index); });
  cls.def("treeIndex", [](const Tree& self, unsigned int index) { return self.treeIndex(); });
}

template<typename Tree, std::size_t i>
auto childAccessor()
{
  return [](Tree& tree) {
    return pybind11::cast(tree.template child(index_constant<i>{}));
  };
}

template<typename Tree, std::size_t... I>
std::array< std::function<pybind11::object(Tree& tree)>, sizeof...(I) >
childAccessors(std::index_sequence<I...>)
{
  static_assert(sizeof...(I) == Tree::degree(),"size of index sequence does not match degree of tree");
  return { (childAccessor<Tree, I>())... };
}

template<typename Tree>
void registerTreeChildAccessor(pybind11::class_<Tree, std::shared_ptr<Tree>>& cls)
{
  const auto accessors = childAccessors<Tree>(std::make_index_sequence<Tree::degree()>{});
  cls.def(
    "child",
    [accessors](Tree& tree, std::size_t i) { return accessors.at(i)(tree); },
    pybind11::arg("i"));
}

template<typename Tree, std::enable_if_t< Tree::isComposite, int > >
void registerTree_(pybind11::handle scope, const char* name)
{
  if( !pybind11::already_registered< Tree >() )
  {
    pybind11::class_<Tree, std::shared_ptr<Tree>> cls(scope, name);
    registerTreeCommon(cls);

    // static variable cls is captured automatically - explicit capture fails with clang:
    //    fatal error: 'cls' cannot be captured because it does not have automatic storage duration
    Hybrid::forEach(std::make_index_sequence<Tree::degree()>{}, [&cls](auto i) {
        using SubTree = typename Tree::template Child<i>::Type;
        std::string subName = std::string("Tree") + std::to_string(i);
        if( !pybind11::already_registered< SubTree >() )
          registerTree_<SubTree>(cls, subName.c_str());
      });
    registerTreeChildAccessor(cls);
  }
}

template<typename Tree, typename = void_t<>>
struct hasFiniteElement
  : std::false_type
{};

template<typename Tree>
struct hasFiniteElement<Tree, void_t<typename Tree::FiniteElement>>
  : std::true_type
{};

template< typename Tree, std::enable_if_t< !hasFiniteElement<Tree>::value, int > = 0>
void registerFiniteElementProperty(pybind11::class_< Tree, std::shared_ptr<Tree> >&)
{
  /* Nothing. */
}

template< typename Tree, std::enable_if_t< hasFiniteElement<Tree>::value, int > = 0>
void registerFiniteElementProperty(pybind11::class_< Tree, std::shared_ptr<Tree> >& cls)
{
  // this should probably be fixed in dune-localfunctions
  if (!pybind11::already_registered< typename Tree::FiniteElement >())
    registerLocalFiniteElement<typename Tree::FiniteElement>(cls);
  cls.def(
    "finiteElement",
    &Tree::finiteElement,
    pybind11::return_value_policy::reference_internal
    );
}

// fatal error: template parameter redefines default argument (clang)
// using static to avoid double registration doesn't work with clang
// because statics remain local to module and are not made unique between modules
// If the class is needed use the TypeRegistry to make this work smoothly
template<typename Tree, std::enable_if_t< Tree::isLeaf, int >>
void registerTree_(pybind11::handle scope, const char* name)
{
  if( !pybind11::already_registered< Tree >() )
  {
    pybind11::class_< Tree, std::shared_ptr<Tree> > cls(scope, name);
    registerTreeCommon(cls);
    registerFiniteElementProperty(cls);
  }
}

template<typename Tree, std::enable_if_t< Tree::isPower, int >>
void registerTree_(pybind11::handle scope, const char* name)
{
  if( !pybind11::already_registered< Tree >() )
  {
    pybind11::class_< Tree, std::shared_ptr<Tree> > cls(scope, name);
    registerTreeCommon(cls);
    registerTree_<typename Tree::ChildType>(cls);
    registerTreeChildAccessor(cls);
  }
}

} /* namespace detail */

template<typename Tree>
void registerTree(pybind11::handle scope, const char* name = "Tree")
{
  detail::registerTree_<Tree>(scope, name);
}

} /* namespace Functions */
} /* namespace Python */
} /* namespace Dune */

#endif
