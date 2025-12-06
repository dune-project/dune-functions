# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

set(HEADERS
  discretefunction.hh
  foreachboundarydof.hh
  globalbasis.hh
  hierarchicvectorwrapper.hh
  interpolate.hh
  subspacebasis.hh
  tree.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/python/functions)
