# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

add_custom_target(build_examples)
add_dependencies(build_tests build_examples)

add_executable("advection-reaction-dg" advection-reaction-dg.cc)
target_link_libraries(advection-reaction-dg PRIVATE Dune::Functions)
set_property(TARGET advection-reaction-dg PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples advection-reaction-dg)

add_executable("interpolation" interpolation.cc)
target_link_libraries(interpolation PRIVATE Dune::Functions)
set_property(TARGET interpolation PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples interpolation)

add_executable("poisson-pq2" poisson-pq2.cc)
target_link_libraries(poisson-pq2 PRIVATE Dune::Functions)
set_property(TARGET poisson-pq2 PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples poisson-pq2)

add_executable("stokes-taylorhood" stokes-taylorhood.cc)
target_link_libraries(stokes-taylorhood PRIVATE Dune::Functions)
set_property(TARGET stokes-taylorhood PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples stokes-taylorhood)

add_executable("poisson-mfem" poisson-mfem.cc)
target_link_libraries(poisson-mfem PRIVATE Dune::Functions)
set_property(TARGET poisson-mfem PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples poisson-mfem)
