# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

add_subdirectory(dune)

dune_python_configure_bindings(
  PATH "."
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
