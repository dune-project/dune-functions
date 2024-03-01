import numpy as np

import dune.grid
from dune.grid import cartesianDomain
from dune.grid import yaspGrid
import dune.functions

# This check will not work currently, because the interface of
# DiscreteGlobalBasisFunction is not provided in python so far.
def checkBasisFunctionInterpolation(basis):
    coeffTol = 1e-10;
    ei = np.zeros(len(basis))
    y = np.zeros(len(basis))
    for i in range(len(basis)):
      ei[i] = 1
      basis.interpolate(y,f)
      if (np.linalg.norm(y-ei) < coeffTol):
          raise ValueError("FE basis function is not interpolated correctly.")
      ei[i] = 0

# Test interpolation of constant function
# Since we cannot deduce a range type automatically,
# a prototype needs to be passed.
def checkConstantInterpolation(basis, rangePrototype):
    coeffTol = 1e-10;
    one = rangePrototype
    one *= 0
    one += 1
    f = lambda x : one
    x = np.zeros(len(basis))
    basis.interpolate(x,f)
    if (np.linalg.norm(x-1) > coeffTol):
        raise ValueError("Interpolation of constant function yields wrong result.")
