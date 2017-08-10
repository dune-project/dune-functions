// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACETYPES_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACETYPES_HH

namespace Dune {
namespace Functions {

  // Define function space types in order to decide what kind of continuity preserving map
  // is used to construct basis shape function on arbitrary elements from a reference element.
  enum FunctionSpaceType {H=0, Hdiv=1, Hcurl=2};

}
}

#endif
