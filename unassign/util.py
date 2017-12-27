from collections import OrderedDict


def uniq(xs):
    """Remove duplicate entries from a list.

    Preserves the order of the input list.
    """
    return list(OrderedDict.fromkeys(xs).keys())

def count_while_equal(xs, val):
    """Count items at start of sequence equal to specified value."""
    n = 0
    for x in xs:
        if x == val:
            n += 1
        else:
            return n
    return n

def count_matching_pairs(xs, gap_char = "-"):
    """Count number of pairs with matching values."""
    match = 0
    total_x1 = 0
    total_x2 = 0
    for x1, x2 in xs:
        if x1 == x2:
            match += 1
        if x1 != gap_char:
            total_x1 += 1
        if x2 != gap_char:
            total_x2 += 1
    return match, total_x1, total_x2
