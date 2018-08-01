# tests that should build and run successfully

# Path to the example grid files in dune-grid
add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${DUNE_GRID_EXAMPLE_GRIDS_PATH}\")

dune_add_test(SOURCES brezzidouglasmarinibasistest.cc)

dune_add_test(SOURCES bsplinebasistest.cc)

dune_add_test(SOURCES gridviewfunctionspacebasistest.cc)

dune_add_test(SOURCES lagrangebasistest.cc)

dune_add_test(SOURCES lagrangedgbasistest.cc)

dune_add_test(SOURCES taylorhoodbasistest.cc)

dune_add_test(SOURCES rannacherturekbasistest.cc)

dune_add_test(SOURCES raviartthomasbasistest.cc)

dune_add_test(SOURCES hierarchicvectorwrappertest.cc)

dune_add_test(SOURCES compositebasistest.cc)

dune_add_test(SOURCES makebasistest.cc)

install(
  FILES
  basistest.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases/test)
