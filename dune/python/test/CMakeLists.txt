# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

dune_python_add_test(NAME pyfunction
                     SCRIPT function.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python)
dune_python_add_test(NAME pydiscretefunction
                     SCRIPT discretefunction.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python)
dune_python_add_test(NAME pycompositebasistest
                     SCRIPT compositebasistest.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python)
dune_python_add_test(NAME pylagrangebasistest
                     SCRIPT lagrangebasistest.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python)
dune_python_add_test(NAME pynedelecbasistest
                     SCRIPT nedelecbasistest.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python)
dune_python_add_test(NAME pyraviartthomasbasistest
                     SCRIPT raviartthomasbasistest.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python)
dune_python_add_test(NAME pysubspacebasistest
                     SCRIPT subspacebasistest.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python)
