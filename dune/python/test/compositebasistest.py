# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

import dune.grid
import dune.functions as functions

import basistest   # The general test suite for function space bases


# Test all currently supported global bases in dimension `dimension`.
# Only the size of the global bases are tested.
def test(dimension):
    lowerLeft = [0] * dimension
    upperRight = [1] * dimension
    elements = [4] * dimension

    grid = dune.grid.structuredGrid(lowerLeft,upperRight,elements)

    # Test a basis consisting of three identical scalar-valued children
    basis = functions.defaultGlobalBasis(grid,
                                         functions.Composite(functions.Lagrange(order=1),
                                                             functions.Lagrange(order=1),
                                                             functions.Lagrange(order=1)))
    basistest.checkBasis(basis)

    # Test a composite of a power node and a Lagrange basis
    basis = functions.defaultGlobalBasis(grid,
                                         functions.Composite(functions.Power(functions.Lagrange(order=2),exponent=dimension),
                                                             functions.Lagrange(order=1)))
    basistest.checkBasis(basis)


# Run tests for grids of dimension 2 and 3
test(2)
test(3)
