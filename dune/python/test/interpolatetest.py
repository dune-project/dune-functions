import numpy as np

import dune.grid
from dune.grid import cartesianDomain
from dune.grid import yaspGrid
import dune.functions

from dune.generator import algorithm
from io import StringIO

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


# Check if indivudual basis functions are interpolated correctly
def checkBasisFunctionInterpolation(basis):
    coeffTol = 1e-10;
    ei = np.zeros(len(basis))
    y = np.zeros(len(basis))
    for i in range(len(basis)):
      ei[i] = 1
      f = asStdFunction(basis.asFunction(ei))
      basis.interpolate(y,f)
      if (np.linalg.norm(y-ei) > coeffTol):
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
