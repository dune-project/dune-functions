# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

import dune.common


# Mark all degrees of freedom on the grid boundary
def forEachBoundaryDOF(basis, callBack):
    """Loop over all degrees of freedom (DOFs) associated to the boundary

       This loops over all DOFs of a basis associated to sub-entities
       on the boundary. This overload will pass a single argument to the
       given loop callback: The global (multi-)index of the boundary DOF.

       Notice that the same DOF may be visited multiple times.

       Parameters
       ----------
       arg0: A function space basis that defines a set of degrees of freedom
       arg1: A callback that will be called with the global index of the visited boundary DOF
    """
    from io import StringIO
    from dune import generator
    code="""
    #include<dune/functions/functionspacebases/boundarydofs.hh>
    template<class Basis, class CallBack>
    void run(const Basis& basis, CallBack& callBack)
    {
      Dune::Functions::forEachBoundaryDOF<Basis,std::function<void(std::size_t)> >(basis, callBack);
    }
    """
    dune.generator.algorithm.run("run",StringIO(code), basis, callBack)
