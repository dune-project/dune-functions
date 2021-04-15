add_subdirectory(functions)

if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory("python")
endif()
