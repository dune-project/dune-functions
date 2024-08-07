# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

dune_add_latex_document(
  SOURCE dune-functions-bases.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR})
dune_add_latex_document(
  SOURCE dune-functions-functions.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR})
