# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

import numpy as np

import dune.grid
from dune.grid import cartesianDomain
from dune.grid import yaspGrid
import dune.functions

from dune.generator import algorithm
from io import StringIO


# Compute a mask vector identifying only the DOFs
# contained in the subspace.
def subspaceMask(basis):
    mask = np.zeros(len(basis))
    localView = basis.localView()
    for element in basis.gridView.elements:
        localView.bind(element)
        tree = localView.tree()
        treeSize = len(tree)
        treeOffset = tree.localIndex(0)
        for i in range(treeOffset, treeOffset + treeSize):
            mask[localView.index(i)[0]] = 1
    return mask

# Manually convert function to std::function
def asStdFunction(f):
    code="""
    #include<utility>
    #include<functional>
    #include<dune/common/fvector.hh>
    template <class F>
    auto run(F f) {
      using Range = typename F::Range;
      using Domain = typename F::Domain;
      return std::function<Range(Domain)>(std::move(f));
    }
    """
    return dune.generator.algorithm.run("run",StringIO(code), f)


# Check if individual basis functions are interpolated correctly
def checkBasisFunctionInterpolation(basis):
    mask = subspaceMask(basis)
    coeffTol = 1e-10;
    ei = np.zeros(len(basis))
    y = np.zeros(len(basis))
    for i in range(len(basis)):
      ei[i] = 1
      ei *= mask
      f = asStdFunction(basis.asFunction(ei))
      basis.interpolate(y,f)
      if (np.linalg.norm(y-ei) > coeffTol):
          raise ValueError("FE basis function is not interpolated correctly.")
      ei[i] = 0


# Test interpolation of constant function
# Since we cannot deduce a range type automatically,
# a prototype needs to be passed.
def checkConstantInterpolation(basis, rangePrototype):
    mask = subspaceMask(basis)
    coeffTol = 1e-10;
    one = rangePrototype
    one *= 0
    one += 1
    f = lambda x : one
    x = np.zeros(len(basis))
    basis.interpolate(x,f)
    if (np.linalg.norm(x-mask) > coeffTol):
        raise ValueError("Interpolation of constant function yields wrong result.")
