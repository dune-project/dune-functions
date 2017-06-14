add_executable("advection-reaction-dg" advection-reaction-dg.cc)
target_link_dune_default_libraries("advection-reaction-dg")

add_executable("poisson-pq2" poisson-pq2.cc)
target_link_dune_default_libraries("poisson-pq2")

add_executable("stokes-taylorhood" stokes-taylorhood.cc)
target_link_dune_default_libraries("stokes-taylorhood")

add_executable("poisson-mfem" poisson-mfem.cc)
target_link_dune_default_libraries("poisson-mfem")

add_executable("poisson-rt0-p0-unstructured" poisson-rt0-p0-unstructured.cc)
target_link_dune_default_libraries("poisson-rt0-p0-unstructured")

add_executable("rt0-basistest" rt0-basistest.cc)
target_link_dune_default_libraries("rt0-basistest")

add_executable("derivative-test" derivative-test.cc)
target_link_dune_default_libraries("derivative-test")
