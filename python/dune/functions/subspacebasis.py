# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt


# A subspace basis for some descendent node, can be obtained
# by passing the indices of the corresponding tree path
# as individual arguments. Passing a tree path at once is
# currently not supported.
#
# Netsed calls to subspaceBasis(subspaceBasis(basis, i), j)
# will create a nested subspace basis. This works fine, but
# is slightly inconsistent with C++ where this is automatically
# resolved to subspaceBasis(basis, i, j).
def subspaceBasis(basis, *args):
    includes = []
    includes += list(basis.cppIncludes)
    includes += ["dune/functions/functionspacebases/subspacebasis.hh"]
    includes += ["dune/python/functions/subspacebasis.hh"]

    # In contrast to C++, the generated wrapper of SubspaceBasis
    # is only constructed from the root basis and no tree path
    # is passed. This is possible, because we only use fully
    # static tree path types that can be default constructed.
    prefixPathTypeName = "Dune::TypeTree::HybridTreePath<"
    prefixPathTypeName += ",".join("Dune::index_constant<"+str(i)+">" for i in args)
    prefixPathTypeName += ">"

    typeName = "Dune::Functions::SubspaceBasis< " + basis.cppTypeName+ "," +  prefixPathTypeName + " >"

    generator = SimpleGenerator("SubspaceBasis", "Dune::Python")

    moduleName = "subspaceBasis_" + hashIt(typeName)
    module = generator.load(includes, typeName, moduleName)
    return module.SubspaceBasis(basis)
