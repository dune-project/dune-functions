# Path to the example grid files in dune-grid
add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${DUNE_GRID_EXAMPLE_GRIDS_PATH}\")

# tests that should build and run successfully
set(TESTS gridviewfunctionspacebasistest taylorhoodbasistest)

set(TESTPROGS ${TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTPROGS})

add_executable("gridviewfunctionspacebasistest"
               gridviewfunctionspacebasistest.cc)
add_executable("taylorhoodbasistest"
               taylorhoodbasistest.cc)
target_link_libraries("gridviewfunctionspacebasistest" "dunecommon" "dunegeometry")
target_link_libraries("taylorhoodbasistest" "dunecommon" "dunegeometry")

add_dune_ug_flags("gridviewfunctionspacebasistest")

# Add the tests to be executed
foreach(_TEST ${TESTPROGS})
  add_test(${_TEST} ${_TEST})
endforeach(_TEST)
