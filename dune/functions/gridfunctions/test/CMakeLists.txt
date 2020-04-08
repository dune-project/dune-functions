# Path to the example grid files in dune-grid
add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${DUNE_GRID_EXAMPLE_GRIDS_PATH}\")

# tests that should build and run successfully
dune_add_test(SOURCES analyticgridviewfunctiontest.cc)

dune_add_test(SOURCES brezzidouglasmarinifunctiontest.cc)

dune_add_test(SOURCES discreteglobalbasisfunctiontest.cc)

dune_add_test(SOURCES gridfunctiontest.cc)

dune_add_test(SOURCES localfunctioncopytest.cc)

dune_add_test(SOURCES raviartthomasfunctiontest.cc)
