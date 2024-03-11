from array import array

# Check consistency of basis.size(prefix)
def checkBasisSizeConsistency(basis, multiIndexSet):

    # Based on the index tree, build a map that contains all possible
    # prefixes and map each prefix to the size (of the subsequent digit).
    prefixSet = dict()

    for index in multiIndexSet:
        prefix = [];
        for digit in index:
            if not tuple(prefix) in prefixSet:
                prefixSet[tuple(prefix)] = 0
            prefixSet[tuple(prefix)] = max(prefixSet[tuple(prefix)], digit+1);
            prefix.append(digit);

        prefixSet[tuple(prefix)] = 0;


    # Now check for all prefixes, if the size computed from the
    # index tree is consistent with basis.size(prefix).
    for prefix, size in prefixSet.items():

        prefixSize = basis.size(prefix);

        if prefixSize != size:
            raise IndexError("basis.size(" + str(prefix) + ") = " + str(prefixSize) + ", but should be " + str(size))


# Check indices of basis:
# - First store the whole index tree in a set
# - Check if this corresponds to a consistent index tree
# - Check if index tree is consistent with basis.size(prefix) and basis.dimension()
def checkBasisIndices(basis):

    multiIndexSet = set()

    localView = basis.localView()
    gridView = basis.gridView

    for element in gridView.elements:
        localView.bind(element)

        for i in range(len(localView)):
            multiIndex = localView.index(i)
            for digit in multiIndex:
                if digit < 0:
                    raise IndexError(str(element) + ": Global multi-index contains negative entry for shape function " + str(i))

            multiIndexSet.add(tuple(multiIndex))

    # TODO: Implement this!
    # checkBasisIndexTreeConsistency(multiIndexSet);
    checkBasisSizeConsistency(basis, multiIndexSet);

    if (basis.dimension != len(multiIndexSet)):
        raise ValueError("basis.dimension() does not equal the total number of basis functions.")


# Check the methods of individual nodes of a function bases tree
# This method calls itself recursively to traverse the entire tree.
def checkTreeNode(node, localIndices):

    if node.isLeaf:

        if node.size() != node.finiteElement.size():
            raise ValueError("Size of leaf node and finite element are different.")

        for i in range(node.size()):

            if node.localIndex(i) >= len(localIndices):
                raise ValueError("Local index exceeds localView size.")

            # Log that we have seen this local index
            localIndices[node.localIndex(i)] += 1

    elif node.isPower or node.isComposite:

        # TODO: Test features of non-leaf nodes

        # Recursively test the children
        for i in range(node.degree()):
            checkTreeNode(node.child(i), localIndices)

    else:
        raise NotImplementedError("Found an unsupported node type")


def checkLocalView(basis, localView, checkLocalIndexSetCompleteness):

    if (localView.size() > localView.maxSize()):
        raise ValueError("localView.size() is " + str(localView.size()) + " but localView.maxSize() is " + str(localView.maxSize()))

    tree1 = localView.tree()
    tree2 = localView.tree()
    if (tree1 != tree2):
        raise ValueError("Multiple calls to localView.tree() do not hand out the same tree.")

    if (tree1.isLeaf and tree1.isPower):
        raise ValueError("Tree claims to be both 'leaf' and 'power'.")

    if (tree1.isLeaf and tree1.isComposite):
        raise ValueError("Tree claims to be both 'leaf' and 'composite'.")

    # Recursively test all nodes in the tree
    localIndices = array('I', [0] * localView.size())

    checkTreeNode(tree1, localIndices)

    # Check that each local index appeared exactly once.
    if checkLocalIndexSetCompleteness:
        for localIndex in localIndices:

            if localIndex < 1:
                raise ValueError("Local index " + str(localIndex) + " did not appear.")
            elif localIndex > 1:
                raise ValueError("Local index " + str(localIndex) + " appeared multiple times.")


# Perform tests that don't modify the basis
def checkConstBasis(basis, *, checkLocalIndexSetCompleteness):

    # Perform all local tests
    localView = basis.localView()
    gridView = basis.gridView

    for element in gridView.elements:
        localView.bind(element)

        # Check the 'element' method
        boundElement = localView.element()
        if boundElement != element:
            raise ValueError("LocalView object does not yield the element it is bound to.")

        checkLocalView(basis,localView, checkLocalIndexSetCompleteness)

    # Perform global index tests
    checkBasisIndices(basis)

    # Perform continuity check
    # TODO


def checkBasis(basis, *, checkLocalIndexSetCompleteness=True):

    # Check the basis
    checkConstBasis(basis, checkLocalIndexSetCompleteness=checkLocalIndexSetCompleteness)

    # Check whether the basis can be copied
    copiedBasis = basis
    checkConstBasis(copiedBasis, checkLocalIndexSetCompleteness=checkLocalIndexSetCompleteness)

    # Can the 'update' method be called?
    gridView = basis.gridView
    # TODO: 'update' method not exposed by the interface yet!
    #basis.update(gridView)
