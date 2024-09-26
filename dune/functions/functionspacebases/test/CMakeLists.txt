# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

# tests that should build and run successfully

# Path to the example grid files in dune-grid
add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${DUNE_GRID_EXAMPLE_GRIDS_PATH}\")

dune_add_test(SOURCES brezzidouglasmarinibasistest.cc LABELS quick)

dune_add_test(SOURCES bsplinebasistest.cc LABELS quick)

dune_add_test(SOURCES containerdescriptortest.cc LABELS quick)

dune_add_test(SOURCES hermitebasistest.cc LABELS quick)

dune_add_test(SOURCES hermitetest.cc LABELS quick)

dune_add_test(SOURCES morleytest.cc LABELS quick)

dune_add_test(SOURCES globalvaluedlfetest.cc LABELS quick)
target_compile_definitions(globalvaluedlfetest PRIVATE -DDUNE_DEPRECATED_INTERPOLATE_CHECK=1)

dune_add_test(SOURCES gridviewfunctionspacebasistest.cc LABELS quick)

dune_add_test(SOURCES lagrangebasistest.cc LABELS quick)

dune_add_test(SOURCES lagrangedgbasistest.cc LABELS quick)

dune_add_test(SOURCES lfebasistest.cc LABELS quick)

dune_add_test(SOURCES nedelecbasistest.cc LABELS quick)

dune_add_test(SOURCES periodicbasistest.cc LABELS quick)

dune_add_test(SOURCES taylorhoodbasistest.cc LABELS quick)

dune_add_test(SOURCES rannacherturekbasistest.cc LABELS quick)

dune_add_test(SOURCES raviartthomasbasistest.cc LABELS quick)

dune_add_test(SOURCES hierarchicvectorwrappertest.cc LABELS quick)

dune_add_test(SOURCES subspacebasistest.cc LABELS quick)

dune_add_test(SOURCES compositebasistest.cc LABELS quick)

dune_add_test(SOURCES makebasistest.cc LABELS quick)

dune_add_test(SOURCES hierarchicallagrangebasistest.cc LABELS quick)

dune_add_test(SOURCES refinedlagrangebasistest.cc
              CMAKE_GUARD Alberta_FOUND)
add_dune_alberta_flags(refinedlagrangebasistest WORLDDIM 2)

install(
  FILES
  basistest.hh
  interpolatetest.hh
  enabledifferentiabilitycheck.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases/test)
