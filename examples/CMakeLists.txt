add_custom_target(build_examples)
add_dependencies(build_tests build_examples)

add_executable("advection-reaction-dg" advection-reaction-dg.cc)
dune_target_enable_all_packages("advection-reaction-dg")
set_property(TARGET advection-reaction-dg PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples advection-reaction-dg)

add_executable("interpolation" interpolation.cc)
dune_target_enable_all_packages("interpolation")
set_property(TARGET interpolation PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples interpolation)

add_executable("poisson-pq2" poisson-pq2.cc)
dune_target_enable_all_packages(poisson-pq2)
set_property(TARGET poisson-pq2 PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples poisson-pq2)

add_executable("stokes-taylorhood" stokes-taylorhood.cc)
dune_target_enable_all_packages("stokes-taylorhood")
set_property(TARGET stokes-taylorhood PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples stokes-taylorhood)

add_executable("poisson-mfem" poisson-mfem.cc)
dune_target_enable_all_packages("poisson-mfem")
set_property(TARGET poisson-mfem PROPERTY EXCLUDE_FROM_ALL 1)
add_dependencies(build_examples poisson-mfem)
