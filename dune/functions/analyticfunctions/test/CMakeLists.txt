# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

# tests that should build and run successfully

dune_add_test(SOURCES monomialsettest.cc LABELS quick)
dune_add_test(SOURCES polynomialtest.cc LABELS quick)
