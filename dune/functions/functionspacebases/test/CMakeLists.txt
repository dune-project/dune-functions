# tests that should build and run successfully
set(TESTS gridviewfunctionspacebasistest)

set(TESTPROGS ${TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTPROGS})

add_executable("gridviewfunctionspacebasistest"
               gridviewfunctionspacebasistest.cc)
target_link_libraries("gridviewfunctionspacebasistest" "dunecommon" "dunegeometry")

# Add the tests to be executed
foreach(_TEST ${TESTPROGS})
  add_test(${_TEST} ${_TEST})
endforeach(_TEST)
