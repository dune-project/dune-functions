# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

import warnings

class Tree(object):
    def __init__(self, name, children=None):
        self.name = name
        self.children = []
        if children is not None:
            assert(all(isinstance(c, Tree) for c in children))
            self.children = list(children)

    def __mul__(self, other):
        return Composite(self, other)

    def __pow__(self, p):
        return Power(self, p)


class Lagrange(Tree):
    def __init__(self, order):
        Tree.__init__(self, "Lagrange")
        self.order = order

    def __repr__(self):
        return "Lagrange<" + str(self.order) + ">"


class DG(Tree):
    def __init__(self, order, dimRange=1):
        Tree.__init__(self, "DG")
        self.order = order
        self.dimRange = dimRange

    def __repr__(self):
        if self.dimRange == 1:
            return "DG<" + str(self.order) + ">"
        else:
            return "DG<" + str(self.order) + ">^" + str(self.dimRange)


class Nedelec(Tree):
    def __init__(self, kind, order):
        Tree.__init__(self, "Nedelec")
        self.kind = kind
        self.order = order

    def __repr__(self):
        return "Nedelec<" + str(self.kind) + "," + str(self.order) + ">"


class RaviartThomas(Tree):
    def __init__(self, order):
        Tree.__init__(self, "RaviartThomas")
        self.order = order

    def __repr__(self):
        return "RaviartThomas<" + str(self.order) + ">"


class Composite(Tree):
    def __init__(self, *args, **kwargs):
        assert len(args) > 0
        Tree.__init__(self, "Composite", args)
        # Please do not change the defaults here without taking the C++ interface into account
        self.blocked = kwargs.get("blocked", False)
        self.layout = kwargs.get("layout", "lexicographic")

    def __repr__(self):
        return "(" + " * ".join(repr(c) for c in self.children) + ")"


class Power(Tree):
    def __init__(self, children, exponent, **kwargs):
        assert children is not None
        Tree.__init__(self, "Power", [children])
        assert len(self.children) == 1
        self.exponent = exponent
        # Please do not change the defaults here without taking the C++ interface into account
        self.blocked = kwargs.get("blocked", False)
        if (not "layout" in kwargs):
            warnings.warn('''The default layout of Power nodes will change from 'lexicographic' to 'interleaved' in release 2.11. Please add the argument 'layout="lexicographic"' explicitly to retain the old behavior!''')
        self.layout = kwargs.get("layout", "lexicographic")

    def __repr__(self):
        if self.exponent == 1:
            return repr(self.children[0])
        else:
            return "[" + repr(self.children[0]) + "]^" + str(self.exponent)
