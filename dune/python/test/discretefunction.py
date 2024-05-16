# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

from dune.grid import structuredGrid
import dune.functions
import numpy as np

# create a simple grid
grid = structuredGrid([0,0],[1,1],[10,10])

# simple P1 basis
basis = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=1))

# boring coefficient vector
coeff = np.zeros( len(basis) )

# discrete function over the grid
f = basis.asFunction(coeff)

# extract the original grid
grid2 = f.grid

# deduce type of discrete function
DiscreteFunction = type(f)

# check constructor of DiscreteFunction
f2 = DiscreteFunction(basis, coeff)

# check if there are still 100 squares
assert( grid2.size(0) == 100 )
