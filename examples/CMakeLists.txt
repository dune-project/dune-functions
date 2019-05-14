add_executable("advection-reaction-dg" advection-reaction-dg.cc)
dune_target_enable_all_packages("advection-reaction-dg")

add_executable("interpolation" interpolation.cc)
dune_target_enable_all_packages("interpolation")

add_executable("poisson-pq2" poisson-pq2.cc)
dune_target_enable_all_packages(poisson-pq2)

add_executable("stokes-taylorhood" stokes-taylorhood.cc)
dune_target_enable_all_packages("stokes-taylorhood")

add_executable("poisson-mfem" poisson-mfem.cc)
dune_target_enable_all_packages("poisson-mfem")
