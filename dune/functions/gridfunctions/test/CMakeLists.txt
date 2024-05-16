# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

# tests that should build and run successfully

dune_add_test(SOURCES analyticgridviewfunctiontest.cc LABELS quick)

dune_add_test(SOURCES composedgridfunctiontest.cc LABELS quick)

dune_add_test(SOURCES discreteglobalbasisfunctiontest.cc LABELS quick)

dune_add_test(SOURCES discreteglobalbasisfunctionderivativetest.cc LABELS quick)

dune_add_test(SOURCES facenormalgridfunctiontest.cc LABELS quick)

dune_add_test(SOURCES gridfunctiontest.cc LABELS quick)

dune_add_test(SOURCES localfunctioncopytest.cc LABELS quick)
