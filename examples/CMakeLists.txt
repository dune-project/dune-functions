add_executable("poisson-pq2" poisson-pq2.cc)
target_link_dune_default_libraries("poisson-pq2")

add_executable("stokes-taylorhood" stokes-taylorhood.cc)
target_link_dune_default_libraries("stokes-taylorhood")
