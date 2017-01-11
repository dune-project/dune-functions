# tests that should build and run successfully

# Path to the example grid files in dune-grid
add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${DUNE_GRID_EXAMPLE_GRIDS_PATH}\")

dune_add_test(SOURCES gridviewfunctionspacebasistest.cc)

dune_add_test(SOURCES taylorhoodbasistest.cc)

dune_add_test(SOURCES hierarchicvectorwrappertest.cc)

dune_add_test(SOURCES compositebasistest.cc)
