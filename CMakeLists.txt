# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.16)
project(dune-functions CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# deactivate global include-directories
dune_policy(SET DP_DEFAULT_INCLUDE_DIRS NEW)

# deactivate global calls to add_dune_all_flags in tests
dune_policy(SET DP_TEST_ADD_ALL_FLAGS NEW)

# start a dune project with information from dune.module
dune_project()
dune_add_library(dunefunctions INTERFACE
  EXPORT_NAME Functions
  LINK_LIBRARIES ${DUNE_LIBS})

if(dune-typetree_FOUND)
  target_link_libraries(dunefunctions INTERFACE Dune::TypeTree)
  dune_mark_module_as_required_dependency(MODULE dune-typetree)
endif()

# set include directories for dunefunctions library
dune_default_include_directories(dunefunctions INTERFACE)

add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("examples")

if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory(python)
endif()


# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project()
