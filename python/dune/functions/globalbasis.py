from __future__ import absolute_import, division, print_function, unicode_literals

from .tree import Composite, DG, Lagrange, Power, Tree

duneFunctionsLayouts = {"lexicographic": "Lexicographic", "interleaved": "Interleaved"}

def basisBuilder(tree):
    assert isinstance(tree, Tree)
    if isinstance(tree, Lagrange):
        scalarBasisBuilder = "Dune::Functions::BasisBuilder::lagrange< " + str(tree.order) + " >()"
        if tree.dimRange != 1:
            return "Dune::Functions::BasisBuilder::power< " + str(tree.dimRange) + " >( " + scalarBasisBuilder + ", Dune::Functions::BasisBuilder::flatInterleaved() )"
        else:
            return scalarBasisBuilder
    elif isinstance(tree, DG):
        raise Exception(repr(tree) + " not supported by dune-functions.")
    elif isinstance(tree, Composite):
        layout = "Dune::Functions::BasisBuilder::" + ("blocked" if tree.blocked else "flat") + duneFunctionsLayouts[tree.layout] + "()"
        return "Dune::Functions::BasisBuilder::composite( " + ", ".join(basisBuilder(c) for c in tree.children) + ", " + layout + " )"
    elif isinstance(tree, Power):
        layout = "Dune::Functions::BasisBuilder::" + ("blocked" if tree.blocked else "flat") + duneFunctionsLayouts[tree.layout] + "()"
        return "Dune::Functions::BasisBuilder::power< " + str(tree.exponent) + " >( " + basisBuilder(tree.children[0]) + ", " + layout + " )"
    else:
        raise Exception("Unknown type of tree: " + repr(tree))


def defaultGlobalBasis(gridView, tree):
    from dune.functions import load

    headers = ["powerbasis", "compositebasis", "lagrangebasis", "subspacebasis", "defaultglobalbasis"]

    includes = []
    includes += list(gridView._includes)
    includes += ["dune/functions/functionspacebases/" + h + ".hh" for h in headers]

    FactoryTag = "decltype( " + basisBuilder(tree) + " )"
    typeName = "Dune::Python::DefaultGlobalBasis< " + gridView._typeName + ", " + FactoryTag + " >"

    return load(includes, typeName).GlobalBasis(gridView)
