from collections import OrderedDict

def uniq(xs):
    """Remove duplicate entries from a list.

    Preserves the order of the input list.
    """
    return OrderedDict.fromkeys(xs).keys()

def count_while_equal(xs, val):
    """Count items at start of sequence equal to specified value."""
    n = 0
    for x in xs:
        if x == val:
            n += 1
        else:
            return n
    return n
