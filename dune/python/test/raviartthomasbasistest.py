# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

import dune.grid
import dune.functions

import basistest   # The general test suite for function space bases

grid = dune.grid.structuredGrid([0,0],[1,1],[5,5])

basis = dune.functions.defaultGlobalBasis(grid, dune.functions.Nedelec(kind=1,order=1))
basistest.checkBasis(basis)
