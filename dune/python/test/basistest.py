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
                    # TODO: Print an element identifier, as in the C++ test
                    raise IndexError("Global multi-index contains negative entry for shape function " + str(i))

            multiIndexSet.add(tuple(multiIndex))

    # TODO: Implement this!
    # checkBasisIndexTreeConsistency(multiIndexSet);
    checkBasisSizeConsistency(basis, multiIndexSet);

    if (basis.dimension != len(multiIndexSet)):
        raise ValueError("basis.dimension() does not equal the total number of basis functions.")


def checkLocalView(basis, localView):

    if (localView.size() > localView.maxSize()):
        raise ValueError("localView.size() is " + str(localView.size()) + " but localView.maxSize() is " + str(localView.maxSize()))

    tree1 = localView.tree()
    tree2 = localView.tree()
    if (tree1 != tree2):
        raise ValueError("Multiple calls to localView.tree() do not hand out the same tree.")

    # TODO: Implement the rest


# Perform tests that don't modify the basis
def checkConstBasis(basis):

    # Perform all local tests
    localView = basis.localView()
    gridView = basis.gridView

    for element in gridView.elements:
        localView.bind(element)

        # Check the 'element' method
        boundElement = localView.element()
        if boundElement != element:
            raise ValueError("LocalView object does not yield the element it is bound to.")

        checkLocalView(basis,localView)

    # Perform global index tests
    checkBasisIndices(basis)

    # Perform continuity check
    # TODO


def checkBasis(basis):

    # Check the basis
    checkConstBasis(basis)

    # Check whether the basis can be copied
    copiedBasis = basis
    checkConstBasis(copiedBasis)

    # Can the 'update' method be called?
    gridView = basis.gridView
    # TODO: 'update' method not exposed by the interface yet!
    #basis.update(gridView)
