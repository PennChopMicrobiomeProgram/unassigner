from collections import OrderedDict


def uniq(xs):
    """Remove duplicate entries from a list.

    Preserves the order of the input list.
    """
    return list(OrderedDict.fromkeys(xs).keys())

