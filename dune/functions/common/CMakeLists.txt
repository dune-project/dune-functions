add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions)
