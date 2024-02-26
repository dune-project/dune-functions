import dune.grid
import dune.functions

import basistest   # The general test suite for function space bases


# Test all currently supported global bases in dimension `dimension`.
# Only the size of the global bases are tested.
def test(dimension):
    lowerLeft = []
    upperRight = []
    elements = []
    for i in range(dimension):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    grid = dune.grid.structuredGrid(lowerLeft,upperRight,elements)

    # Test a first-order Lagrange basis
    basisLagrange1 = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=1))
    basistest.checkBasis(basisLagrange1)

    # Test a second-order Lagrange basis
    basisLagrange2 = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=2))
    basistest.checkBasis(basisLagrange2)

    # Test a vector-value Lagrange basis
    basisPower = dune.functions.defaultGlobalBasis(grid, dune.functions.Power(dune.functions.Lagrange(order=1),exponent=dimension))
    assert(len(basisPower) == len(basisLagrange1)*dimension)

    # Test a vector-value Lagrange basis with different orders
    basisComposite = dune.functions.defaultGlobalBasis(grid, dune.functions.Composite(dune.functions.Power(dune.functions.Lagrange(order=2),exponent=dimension),
                                                                                      dune.functions.Lagrange(order=1)))
    assert(len(basisComposite) == len(basisLagrange2)*dimension + len(basisLagrange1))

# Run tests for grids of dimension 2 and 3
test(2)
test(3)
