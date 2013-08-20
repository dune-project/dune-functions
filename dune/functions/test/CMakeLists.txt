set(TESTS
    hello)

add_executable(hello hello.cc)
target_link_libraries(hello dunefunctions)

add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})

foreach(_test ${TESTS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})
