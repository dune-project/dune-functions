# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dune.common.checkconfiguration import assertHave, ConfigurationError

from .boundarydofs import *
from .globalbasis import defaultGlobalBasis
from .subspacebasis import subspaceBasis
from .tree import *

registry = dict()
registry["globalBasis"] = {
        "default" : defaultGlobalBasis
    }

generator = SimpleGenerator("GlobalBasis", "Dune::Python")

def load(includes, typeName, *args):
    includes = includes + ["dune/python/functions/globalbasis.hh"]
    moduleName = "globalbasis_" + hashIt(typeName)
    module = generator.load(includes, typeName, moduleName, *args)
    return module
