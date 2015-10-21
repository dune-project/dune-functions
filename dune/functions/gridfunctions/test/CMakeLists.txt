# tests that should build and run successfully
set(TESTS analyticgridviewfunctiontest discreteglobalbasisfunctiontest gridfunctiontest)

set(TESTPROGS ${TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTPROGS})

add_executable("analyticgridviewfunctiontest"
               analyticgridviewfunctiontest.cc)
target_link_libraries("analyticgridviewfunctiontest" ${DUNE_LIBS})

add_executable("discreteglobalbasisfunctiontest"
               discreteglobalbasisfunctiontest.cc)
target_link_libraries("discreteglobalbasisfunctiontest" ${DUNE_LIBS})

add_executable("gridfunctiontest"
               gridfunctiontest.cc)
target_link_libraries("gridfunctiontest" ${DUNE_LIBS})

# Add the tests to be executed
foreach(_TEST ${TESTPROGS})
  add_test(${_TEST} ${_TEST})
endforeach(_TEST)
