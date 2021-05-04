#ifndef DUNE_PYTHON_FUNCTIONS_TREE_HH
#define DUNE_PYTHON_FUNCTIONS_TREE_HH

#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>

// for void_t, in <type_traits> for C++17
#include <dune/common/typetraits.hh>
#include <dune/common/visibility.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/localfunctions/localfiniteelement.hh>

namespace Dune {
namespace Python {
namespace Functions {

namespace detail {

template<typename Tree, std::enable_if_t< Tree::isComposite, int > = 0>
pybind11::class_<Tree> registerTree_(pybind11::handle scope, const char* name = "Tree");

template<typename Tree, std::enable_if_t< Tree::isLeaf, int > = 0>
pybind11::class_<Tree> registerTree_(pybind11::handle scope, const char* name = "Tree");

template<typename Tree, std::enable_if_t< Tree::isPower, int > = 0>
pybind11::class_<Tree> registerTree_(pybind11::handle scope, const char* name = "Tree");

template<typename Tree>
void registerTreeCommon(pybind11::class_<Tree>& cls)
{
  /* dune-typetree properties */
  cls.def_property_readonly_static("isComposite", [](pybind11::object) { return Tree::isComposite; });
  cls.def_property_readonly_static("isLeaf", [](pybind11::object) { return Tree::isLeaf; });
  cls.def_property_readonly_static("isPower", [](pybind11::object) { return Tree::isPower; });
  cls.def_property_readonly_static("degree", [](pybind11::object) { return Tree::degree(); });
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
DUNE_EXPORT void registerTreeChildAccessor(pybind11::class_<Tree>& cls)
{
  static const auto accessors = childAccessors<Tree>(std::make_index_sequence<Tree::degree()>{});
  cls.def(
    "__getitem__",
    [](Tree& tree, std::size_t i) { return accessors.at(i)(tree); },
    pybind11::arg("i"));
}

template<typename Tree, std::enable_if_t< Tree::isComposite, int > >
DUNE_EXPORT pybind11::class_<Tree> registerTree_(pybind11::handle scope, const char* name)
{
  static pybind11::class_<Tree> cls(scope, name);
  registerTreeCommon(cls);

  Hybrid::forEach(std::make_index_sequence<Tree::degree()>{}, [&cls](auto i) {
      using SubTree = typename Tree::template Child<i>::Type;
      std::string subName = std::string("Tree") + std::to_string(i);
      if( !pybind11::already_registered< SubTree >() )
        registerTree_<SubTree>(cls, subName.c_str());
    });
  registerTreeChildAccessor(cls);

  return cls;
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
void registerFiniteElementProperty(pybind11::class_< Tree >&)
{
  /* Nothing. */
}

template< typename Tree, std::enable_if_t< hasFiniteElement<Tree>::value, int > = 0>
void registerFiniteElementProperty(pybind11::class_< Tree >& cls)
{
  registerLocalFiniteElement<typename Tree::FiniteElement>(cls);
  cls.def_property_readonly(
    "finiteElement",
    [](const Tree& tree) { return &tree.finiteElement(); },
    pybind11::return_value_policy::reference_internal
    );
}

template<typename Tree, std::enable_if_t< Tree::isLeaf, int > = 0>
DUNE_EXPORT pybind11::class_<Tree> registerTree_(pybind11::handle scope, const char* name)
{
  static pybind11::class_< Tree > cls(scope, name);
  if( !pybind11::already_registered< Tree >() )
    registerTreeCommon(cls);

  registerFiniteElementProperty(cls);

  // register localIndex
  cls.def("localIndex", [](Tree& tree, unsigned int index) { return tree.localIndex(index); });
  return cls;
}

template<typename Tree, std::enable_if_t< Tree::isPower, int > = 0>
DUNE_EXPORT pybind11::class_<Tree> registerTree_(pybind11::handle scope, const char* name)
{
  static pybind11::class_< Tree > cls(scope, name);
  registerTreeCommon(cls);

  if( !pybind11::already_registered< typename Tree::ChildType >() )
    registerTree_<typename Tree::ChildType>(cls);
  registerTreeChildAccessor(cls);

  return cls;
}

} /* namespace detail */

template<typename Tree>
DUNE_EXPORT pybind11::class_<Tree> registerTree(pybind11::handle scope, const char* name = "Tree")
{
  static auto cls = detail::registerTree_<Tree>(scope, name);
  return cls;
}

} /* namespace Functions */
} /* namespace Python */
} /* namespace Dune */

#endif
