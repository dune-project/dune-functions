# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

import numpy as np

import dune.grid
import dune.functions as functions

import basistest        # The general test suite for function space bases
import interpolatetest  # The general test suite for interpolation into bases


# Test all currently supported global bases in dimension `dimension`.
# Only the size of the global bases are tested.
def test(dimension):
    lowerLeft = [-1] * dimension
    upperRight = [1] * dimension
    elements = [3] * dimension

    grid = dune.grid.structuredGrid(lowerLeft,upperRight,elements)

    # Test a first-order Lagrange basis
    basis1 = functions.defaultGlobalBasis(grid, functions.Lagrange(order=1))
    basistest.checkBasis(basis1)

    # Test a second-order Lagrange basis
    basis2 = functions.defaultGlobalBasis(grid, functions.Lagrange(order=2))
    basistest.checkBasis(basis2)

    # Test a vector-value Lagrange basis
    basisPower = functions.defaultGlobalBasis(grid, functions.Power(functions.Lagrange(order=1),exponent=dimension))
    basistest.checkBasis(basisPower)

    # Test interpolation of constant functions
    interpolatetest.checkConstantInterpolation(basis1, 0)
    interpolatetest.checkConstantInterpolation(basis2, 0)
    interpolatetest.checkConstantInterpolation(basisPower, np.zeros(dimension))

    # Test interpolation of all basis functions
    interpolatetest.checkBasisFunctionInterpolation(basis1)
    interpolatetest.checkBasisFunctionInterpolation(basis2)
    interpolatetest.checkBasisFunctionInterpolation(basisPower)

# Run tests for grids of dimension 2 and 3
test(2)
test(3)
